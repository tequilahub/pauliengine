"""Run the PauliCycle commutator benchmark over powers of two, fit a power law
to each curve, and produce annotated plots — once for dense base strings and
once for bounded-weight ones.

A power law is  T(n) = C * n**alpha ; on log-log axes it is a straight line of
slope alpha. We fit alpha by least squares on log T vs log n, using only the
upper half of the sweep (the asymptotic regime) so the constant per-call
overhead at small n does not bias the exponent.

Run:

    python -m tools.benchmarks.cycle_fit_plot
"""

from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

# Support both invocation styles:
#   python -m tools.benchmarks.cycle_fit_plot   (run as a package module)
#   python cycle_fit_plot.py                    (run as a plain script)
if __package__:
    from . import hardware
    from .cycle_benchmark import measure_point
else:
    # Plain-script mode: put the repo root on sys.path and import via the
    # package, so intra-package imports resolve too.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from tools.benchmarks import hardware
    from tools.benchmarks.cycle_benchmark import measure_point


def run_sweep(exponents: list[int], weight: int | None,
              repeats: int, warmup: int, seed: int) -> dict:
    measurements = []
    for e in exponents:
        n = 2 ** e
        print(f"[fit-bench] weight={weight}  n={n}", flush=True)
        measurements.append(measure_point(n, repeats, warmup, seed, weight))
    return {
        "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "config": {"exponents": exponents, "weight": weight,
                   "repeats": repeats, "warmup": warmup, "seed": seed},
        "hardware": hardware.collect(),
        "measurements": measurements,
    }


def power_law_fit(n: np.ndarray, y: np.ndarray, tail_from: int):
    """Least-squares fit log y = alpha*log n + log C on the points n >= tail_from.

    Returns (alpha, C, r2, mask) where mask selects the fitted points.
    """
    mask = n >= tail_from
    x, yy = np.log(n[mask]), np.log(y[mask])
    alpha, b = np.polyfit(x, yy, 1)
    yhat = alpha * x + b
    ss_res = np.sum((yy - yhat) ** 2)
    ss_tot = np.sum((yy - yy.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return alpha, float(np.exp(b)), r2, mask


def plot_with_fit(payload: dict, out_path: Path) -> dict:
    import matplotlib.pyplot as plt

    ms = sorted(payload["measurements"], key=lambda m: m["n_qubits"])
    n = np.array([m["n_qubits"] for m in ms], float)
    weight = payload["config"]["weight"]
    # Fit from n = 2**5 = 32 onward (skip the small-n, overhead-dominated points).
    tail_from = 32.0

    fig, ax = plt.subplots(figsize=(8, 6))
    fits = {}
    colors = {"naive": "tab:blue", "cycle": "tab:orange"}
    labels = {"naive": "QubitHamiltonian commutator",
              "cycle": "PauliCycleSum"}
    for key in ("naive", "cycle"):
        y = np.array([m[f"{key}_time_mean_s"] for m in ms], float)
        yerr = np.array([m[f"{key}_time_stdev_s"] for m in ms], float)
        alpha, C, r2, mask = power_law_fit(n, y, tail_from)
        fits[key] = {"alpha": alpha, "C": C, "r2": r2}
        ax.errorbar(n, y, yerr=yerr, marker="o", capsize=3,
                    color=colors[key],
                    label=f"{labels[key]}\n   fit  α={alpha:.2f}")
        # Overlay the fitted power law over the fitted range.
        xf = np.array([n[mask].min(), n[mask].max()])
        ax.plot(xf, C * xf ** alpha, ls="--", color=colors[key], alpha=0.7)

    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("n_qubits  (= number of rotations)")
    ax.set_ylabel("commutator time [s] (mean ± stdev)")
    weight_label = "dense" if weight is None else f"weight {weight}"
    ax.set_title(f"Power-law fit — {weight_label}  (fit on n ≥ {int(tail_from)})")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(fontsize=8, loc="upper left")

    # Write the fitted power laws directly onto the plot.
    lines = [r"power-law fit:  $T = C \cdot n^{\alpha}$"]
    for key, name in (("naive", "naive"), ("cycle", "cycle")):
        f = fits[key]
        lines.append(f"{name}:  α = {f['alpha']:.2f},  C = {f['C']:.2e}")
    ax.text(0.97, 0.03, "\n".join(lines), transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="gray"))

    hw_label = hardware.short_label(payload.get("hardware", {}))
    fig.text(0.5, 0.01, hw_label, ha="center", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return fits


def main() -> int:
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    repeats, warmup, seed = 30, 3, 0

    # (name, exponents, weight)
    runs = [
        ("dense", list(range(3, 12)), None),   # n = 8 .. 2048
        ("weight3", list(range(3, 13)), 3),    # n = 8 .. 4096
    ]

    summary = {}
    for name, exps, weight in runs:
        payload = run_sweep(exps, weight, repeats, warmup, seed)
        json_path = results_dir / f"cycle-fit-{name}.json"
        png_path = results_dir / f"cycle-fit-{name}.png"
        json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        fits = plot_with_fit(payload, png_path)
        summary[name] = fits
        print(f"  -> {png_path}")

    print("\n=== Power-law exponents (fit on n >= 32) ===")
    for name, fits in summary.items():
        print(f"  {name}:")
        for key, f in fits.items():
            print(f"      {key:6}  alpha = {f['alpha']:.2f}   (R2 = {f['r2']:.3f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
