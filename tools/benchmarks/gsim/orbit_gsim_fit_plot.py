"""Run the Pauli-orbit commutator benchmark for three routes (naive, PauliEngine
orbit, g-sim), fit a power law to each curve, and produce annotated plots.

A power law is  T(n) = C * n**alpha ; on log-log axes it is a straight line of
slope alpha. We fit alpha by least squares on log T vs log n over the whole sweep
(all points, from the smallest n).

g-sim is used automatically if the ``gsim`` package is installed. If unavailable,
only the naive and orbit curves are produced.

Run:

    python -m tools.benchmarks.orbit_gsim_fit_plot
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

if __package__:
    from ..utils import hardware
    from ..utils.benchmark import has_route, route_stats
    from ..utils.gsim_compat import load_orbits
    from .orbit_gsim_benchmark import measure_point
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from tools.benchmarks.utils import hardware
    from tools.benchmarks.utils.benchmark import has_route, route_stats
    from tools.benchmarks.utils.gsim_compat import load_orbits
    from tools.benchmarks.gsim.orbit_gsim_benchmark import measure_point


def run_sweep(n_values: list[int], weight: int, repeats: int, warmup: int,
              seed: int, gsim) -> dict:
    measurements = []
    for n in n_values:
        print(f"[orbit-gsim-fit] weight={weight}  n={n}", flush=True)
        measurements.append(measure_point(n, repeats, warmup, seed, weight, gsim))
    return {
        "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "config": {"n_values": n_values, "weight": weight, "repeats": repeats,
                   "warmup": warmup, "seed": seed, "gsim": gsim is not None},
        "hardware": hardware.collect(),
        "measurements": measurements,
    }


def power_law_fit(n: np.ndarray, y: np.ndarray, tail_from: float):
    """Least-squares fit log y = alpha*log n + log C on the points n >= tail_from."""
    mask = n >= tail_from
    x, yy = np.log(n[mask]), np.log(y[mask])
    alpha, b = np.polyfit(x, yy, 1)
    yhat = alpha * x + b
    ss_res = np.sum((yy - yhat) ** 2)
    ss_tot = np.sum((yy - yy.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return alpha, float(np.exp(b)), r2, mask


def plot_with_fit(payload: dict, out_path: Path, tail_from: float) -> dict:
    import matplotlib.pyplot as plt

    ms = sorted(payload["measurements"], key=lambda m: m["n_qubits"])
    n = np.array([m["n_qubits"] for m in ms], float)
    weight = payload["config"]["weight"]
    has_gsim = has_route(ms, "gsim")

    routes = [("naive", "QubitHamiltonian commutator", "tab:blue"),
              ("orbit", "PauliEngine PauliOrbit", "tab:orange")]
    if has_gsim:
        routes.append(("gsim", "g-sim orbit_commutator", "tab:green"))

    fig, ax = plt.subplots(figsize=(8, 6))
    fits = {}
    for key, label, color in routes:
        y, yerr = (np.array(v, float) for v in zip(*(route_stats(m, key) for m in ms)))
        alpha, C, r2, mask = power_law_fit(n, y, tail_from)
        fits[key] = {"alpha": alpha, "C": C, "r2": r2}
        ax.errorbar(n, y, yerr=yerr, marker="o", capsize=3, color=color,
                    label=f"{label}\n   fit  α={alpha:.2f}")
        xf = np.array([n[mask].min(), n[mask].max()])
        ax.plot(xf, C * xf ** alpha, ls="--", color=color, alpha=0.7)

    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("n_qubits")
    ax.set_ylabel("commutator time [s] (geo. mean ± geo. spread)")
    fit_range = "over all n" if tail_from <= n.min() else f"on n ≥ {int(tail_from)}"
    ax.set_title(f"Power-law fit — orbit weight {weight}  (fit {fit_range})")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(fontsize=8, loc="upper left")

    lines = [r"power-law fit:  $T = C \cdot n^{\alpha}$"]
    for key, _, _ in routes:
        f = fits[key]
        lines.append(f"{key}:  α = {f['alpha']:.2f},  C = {f['C']:.2e}")
    ax.text(0.97, 0.03, "\n".join(lines), transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="gray"))

    hw_label = hardware.short_label(payload.get("hardware", {}))
    fig.text(0.5, 0.01, hw_label, ha="center", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return fits


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="tools.benchmarks.orbit_gsim_fit_plot",
        description="Power-law fit of the orbit commutator: PauliEngine vs g-sim vs naive.",
    )
    p.add_argument("--repeats", type=int, default=15)
    p.add_argument("--warmup", type=int, default=3)
    p.add_argument("--seed", type=int, default=0)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    gsim = load_orbits()
    if gsim is None:
        print("[orbit-gsim-fit] g-sim not available — fitting naive + orbit only.",
              file=sys.stderr)

    results_dir = Path(__file__).resolve().parents[1] / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # (name, n-values, weight, tail_from). Fit from the first point (tail_from=0);
    # kept small (orbits expand as ~n^weight).
    runs = [
        ("weight2", [4, 6, 8, 12, 16, 24, 32, 48, 64], 2, 0.0),
    ]

    summary = {}
    for name, n_values, weight, tail_from in runs:
        payload = run_sweep(n_values, weight, args.repeats, args.warmup, args.seed, gsim)
        json_path = results_dir / f"orbit-gsim-fit-{name}.json"
        png_path = results_dir / f"orbit-gsim-fit-{name}.png"
        json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        fits = plot_with_fit(payload, png_path, tail_from)
        summary[name] = fits
        print(f"  -> {png_path}")

    print("\n=== Power-law exponents ===")
    for name, fits in summary.items():
        print(f"  {name}:")
        for key, f in fits.items():
            print(f"      {key:6}  alpha = {f['alpha']:.2f}   (R2 = {f['r2']:.3f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
