"""Scaling benchmark: PauliEngine orbit commutator vs g-sim — no naive route.

Unlike ``orbit_gsim_benchmark`` this drops the expanded-QubitHamiltonian route
entirely. Materialising an orbit into explicit Pauli strings is what limits the
qubit count (it grows as ~n^weight), so removing it lets both poly(n) routes —

* **orbit**  — PauliEngine's ``PauliOrbit.commutator`` (C++)
* **gsim**   — g-sim's ``orbit_commutator`` Numba kernel

— be pushed to much larger ``n``. Correctness is still cross-checked once per
point by comparing the set of output orbit labels (normalization independent),
which is cheap and needs no expansion.

g-sim is required here (the whole point is the comparison); if it is not
installed the run aborts. It is loaded automatically via ``utils.gsim_compat``.

Both implementations tabulate factorials in float64, which overflows at 171!, so
``n`` above ~170 is skipped (see ``_MAX_FACTORIAL_N``). That still reaches far
higher qubit counts than the naive route (limited to ~64 by string expansion).

Run (from the repo root):

    python -m tools.benchmarks.gsim.orbit_gsim_scale \
        --n-qubits 8 16 32 64 128 --weight 2 \
        --repeats 10 --warmup 3
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import random
import sys
from pathlib import Path

import numpy as np

if __package__:
    from ..utils import hardware
    from ..utils.benchmark import _time_call, route_stats
    from ..utils.gsim_compat import load_orbits
    from ..native.orbit_benchmark import n_terms, random_label
    from .orbit_gsim_benchmark import _gsim_labels
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from tools.benchmarks.utils import hardware
    from tools.benchmarks.utils.benchmark import _time_call, route_stats
    from tools.benchmarks.utils.gsim_compat import load_orbits
    from tools.benchmarks.native.orbit_benchmark import n_terms, random_label
    from tools.benchmarks.gsim.orbit_gsim_benchmark import _gsim_labels

from pauliengine import PauliOrbit

# Both implementations tabulate factorials as float64, which overflows at 171!
# (> 1.8e308). Beyond this the coefficient arithmetic breaks (g-sim raises, and
# PauliEngine's coefficients go to nan) even though the algorithm itself is
# poly(n). This caps the meaningful sweep — still far beyond the naive route.
_MAX_FACTORIAL_N = 170


def measure_point(n_qubits: int, repeats: int, warmup: int, seed: int,
                  weight: int, gsim) -> dict:
    """Time the orbit and g-sim routes at one qubit count and check labels agree."""
    orbit_commutator, precompute = gsim

    rng = random.Random(seed + n_qubits)
    label_a = random_label(n_qubits, rng, weight)
    label_b = random_label(n_qubits, rng, weight)

    a = PauliOrbit(n_qubits, *label_a)
    b = PauliOrbit(n_qubits, *label_b)
    orbit_call = lambda: a.commutator(b)

    fact = precompute(n_qubits)
    pa, qa, ra = label_a
    pb, qb, rb = label_b
    gsim_call = lambda: orbit_commutator(n_qubits, pa, qa, ra, pb, qb, rb, fact)

    # Label-set correctness (no expansion): both must yield the same output orbits.
    orbit_labels = {(o.p, o.q, o.r) for o in a.commutator(b).data}
    if _gsim_labels(*gsim_call()) != orbit_labels:
        raise AssertionError(f"gsim vs orbit label mismatch at n_qubits={n_qubits}")

    return {
        "n_qubits": n_qubits,
        "label_a": list(label_a),
        "label_b": list(label_b),
        "orbit_strings": n_terms(n_qubits, *label_a),
        "n_output_orbits": len(orbit_labels),
        "repeats": repeats,
        "orbit_times_s": _time_call(orbit_call, repeats=repeats, warmup=warmup),
        "gsim_times_s": _time_call(gsim_call, repeats=repeats, warmup=max(warmup, 2)),
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


def plot(payload: dict, out_path: Path, tail_from: float = 0.0) -> dict:
    import matplotlib.pyplot as plt

    ms = sorted(payload["measurements"], key=lambda m: m["n_qubits"])
    n = np.array([m["n_qubits"] for m in ms], float)
    weight = payload["config"]["weight"]

    fig, ax = plt.subplots(figsize=(8, 6))
    fits = {}
    for key, label, color in (("orbit", "PauliEngine PauliOrbit", "tab:orange"),
                              ("gsim", "g-sim orbit_commutator", "tab:green")):
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
    ax.set_title(f"Orbit commutator scaling — PauliEngine vs g-sim, weight {weight}  (fit {fit_range})")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(fontsize=8, loc="upper left")

    lines = [r"power-law fit:  $T = C \cdot n^{\alpha}$"]
    for key in ("orbit", "gsim"):
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
        prog="tools.benchmarks.gsim.orbit_gsim_scale",
        description="Scaling benchmark of the orbit commutator: PauliEngine vs g-sim (no naive route).",
    )
    p.add_argument("--n-qubits", type=int, nargs="+",
                   default=[8, 16, 24, 32, 48, 64, 96, 128])
    p.add_argument("--weight", type=int, default=2,
                   help="Total non-identity weight p+q+r of each orbit label.")
    p.add_argument("--repeats", type=int, default=10)
    p.add_argument("--warmup", type=int, default=3)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--output", type=Path, default=None)
    p.add_argument("--label", type=str, default=None)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    gsim = load_orbits()
    if gsim is None:
        print("[orbit-gsim-scale] g-sim is required for this benchmark but is not "
              "available — aborting.", file=sys.stderr)
        return 1

    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    output = args.output or (Path(__file__).resolve().parents[1] / "results" / f"orbit-gsim-scale-{stamp}.json")
    output.parent.mkdir(parents=True, exist_ok=True)

    measurements = []
    for n in args.n_qubits:
        if n > _MAX_FACTORIAL_N:
            print(f"[orbit-gsim-scale] skipping n={n}: float64 factorials overflow "
                  f"beyond n={_MAX_FACTORIAL_N}", file=sys.stderr, flush=True)
            continue
        print(f"[orbit-gsim-scale] n_qubits={n} weight={args.weight}", file=sys.stderr, flush=True)
        measurements.append(measure_point(n, args.repeats, args.warmup, args.seed, args.weight, gsim))

    payload = {
        "schema_version": 2,
        "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "label": args.label,
        "benchmark": "pauli_orbit_commutator_scaling_vs_gsim",
        "config": {
            "n_qubits_values": args.n_qubits,
            "weight": args.weight,
            "repeats": args.repeats,
            "warmup": args.warmup,
            "seed": args.seed,
        },
        "hardware": hardware.collect(),
        "measurements": measurements,
    }
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    png_path = output.with_suffix(".png")
    fits = plot(payload, png_path)

    print(f"wrote {len(measurements)} measurements → {output}", file=sys.stderr)
    print(f"wrote plot → {png_path}", file=sys.stderr)
    print("\n=== Power-law exponents ===", file=sys.stderr)
    for key, f in fits.items():
        print(f"  {key:6}  alpha = {f['alpha']:.2f}   (R2 = {f['r2']:.3f})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
