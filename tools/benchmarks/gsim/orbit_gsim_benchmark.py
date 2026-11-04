"""Benchmark: Pauli-orbit commutator — three routes compared.

Same setup as ``orbit_benchmark`` but with g-sim added as a third route:

* **naive**  — expand each orbit into an explicit QubitHamiltonian and commute
  the two full operators (``QubitHamiltonian.commutator``).
* **orbit**  — PauliEngine's Pauli-orbit commutator (``PauliOrbit.commutator``),
  the contingency-table formula in C++.
* **gsim**   — g-sim's ``orbit_commutator`` Numba kernel
  (https://github.com/adelina-b/g-sim, src/gsim/gorbits.py).

All three enumerate the same contingency tables and therefore produce the same
set of output orbit labels (checked once per point; g-sim returns raw
combinatorial coefficients without PauliEngine's 2/(N1 N2) prefactor, so only the
labels are compared, which is normalization independent).

g-sim is optional and used automatically if the ``gsim`` package is installed
(see ``gsim_compat`` for the PauliEngine import shim). Otherwise the g-sim route
is skipped.

Usage (from the repo root):

    python -m tools.benchmarks.orbit_gsim_benchmark \
        --n-qubits 4 6 8 12 16 24 32 --weight 2 \
        --repeats 5 --warmup 2 \
        --output tools/benchmarks/results/orbit-gsim-latest.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import random
import sys
from pathlib import Path

# Support both invocation styles.
if __package__:
    from ..utils import hardware
    from ..utils.benchmark import _time_call, has_route, route_stats
    from ..utils.gsim_compat import load_orbits
    from ..native.orbit_benchmark import (
        n_terms,
        orbit_to_hamiltonian,
        random_label,
        _orbitsum_to_hamiltonian,
        _qh_close,
    )
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from tools.benchmarks.utils import hardware
    from tools.benchmarks.utils.benchmark import _time_call, has_route, route_stats
    from tools.benchmarks.utils.gsim_compat import load_orbits
    from tools.benchmarks.native.orbit_benchmark import (
        n_terms,
        orbit_to_hamiltonian,
        random_label,
        _orbitsum_to_hamiltonian,
        _qh_close,
    )

from pauliengine import PauliOrbit


def _gsim_labels(kp, kq, kr, kv, tol: float = 1e-9) -> set:
    """Set of output labels with non-zero coefficient from a g-sim result."""
    return {
        (int(kp[i]), int(kq[i]), int(kr[i]))
        for i in range(len(kv))
        if abs(kv[i]) > tol
    }


def measure_point(n_qubits: int, repeats: int, warmup: int, seed: int,
                  weight: int, gsim=None) -> dict:
    """Time the routes at one qubit count and assert they agree (labels)."""
    rng = random.Random(seed + n_qubits)
    label_a = random_label(n_qubits, rng, weight)
    label_b = random_label(n_qubits, rng, weight)

    a = PauliOrbit(n_qubits, *label_a)
    b = PauliOrbit(n_qubits, *label_b)

    h_a = orbit_to_hamiltonian(n_qubits, *label_a)
    h_b = orbit_to_hamiltonian(n_qubits, *label_b)
    naive_call = lambda: h_a.commutator(h_b)
    orbit_call = lambda: a.commutator(b)

    # Correctness: orbit route vs the expanded operator.
    if not _qh_close(_orbitsum_to_hamiltonian(a.commutator(b)), naive_call()):
        raise AssertionError(f"orbit vs naive mismatch at n_qubits={n_qubits}")

    record = {
        "n_qubits": n_qubits,
        "label_a": list(label_a),
        "label_b": list(label_b),
        "orbit_strings": n_terms(n_qubits, *label_a),
        "result_terms": len(naive_call()),
        "repeats": repeats,
    }

    for name, call in (("naive", naive_call), ("orbit", orbit_call)):
        record[f"{name}_times_s"] = _time_call(call, repeats=repeats, warmup=warmup)

    if gsim is not None:
        orbit_commutator, precompute = gsim
        fact = precompute(n_qubits)
        pa, qa, ra = label_a
        pb, qb, rb = label_b
        gsim_call = lambda: orbit_commutator(n_qubits, pa, qa, ra, pb, qb, rb, fact)

        # Label-set correctness check (normalization independent).
        orbit_labels = {(o.p, o.q, o.r) for o in a.commutator(b).data}
        if _gsim_labels(*gsim_call()) != orbit_labels:
            raise AssertionError(f"gsim vs orbit label mismatch at n_qubits={n_qubits}")

        # extra warmup so the Numba JIT compile is never inside the timed region.
        record["gsim_times_s"] = _time_call(gsim_call, repeats=repeats, warmup=max(warmup, 2))

    return record


def plot(payload: dict, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    ms = sorted(payload["measurements"], key=lambda m: m["n_qubits"])
    xs = [m["n_qubits"] for m in ms]
    weight = payload.get("config", {}).get("weight")
    has_gsim = has_route(ms, "gsim")

    fig, ax = plt.subplots(figsize=(8, 6))
    routes = [
        ("naive", "naive: full QubitHamiltonian commutator", "s", "tab:blue"),
        ("orbit", "orbit: PauliEngine PauliOrbit", "o", "tab:orange"),
    ]
    if has_gsim:
        routes.append(("gsim", "gsim: g-sim orbit_commutator (numba)", "^", "tab:green"))

    for key, label, marker, color in routes:
        ys, yerr = zip(*(route_stats(m, key) for m in ms))
        ax.errorbar(xs, ys, yerr=yerr, marker=marker, capsize=3, color=color, label=label)

    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("n_qubits")
    ax.set_ylabel("commutator time [s] (geo. mean ± geo. spread)")
    ax.set_title(f"Pauli-orbit commutator: PauliEngine vs g-sim vs naive (weight {weight})")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(fontsize=8)

    hw_label = hardware.short_label(payload.get("hardware", {}))
    fig.text(0.5, 0.01, hw_label, ha="center", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="tools.benchmarks.orbit_gsim_benchmark",
        description="Benchmark the Pauli-orbit commutator: PauliEngine vs g-sim vs naive.",
    )
    p.add_argument("--n-qubits", type=int, nargs="+", default=[4, 6, 8, 12, 16, 24, 32])
    p.add_argument("--weight", type=int, default=2,
                   help="Total non-identity weight p+q+r of each orbit label.")
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--warmup", type=int, default=2)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--output", type=Path, default=None)
    p.add_argument("--label", type=str, default=None)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    gsim = load_orbits()
    if gsim is None:
        print("[orbit-gsim] g-sim not available — measuring naive + orbit only.",
              file=sys.stderr)

    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    output = args.output or (Path(__file__).resolve().parents[1] / "results" / f"orbit-gsim-{stamp}.json")
    output.parent.mkdir(parents=True, exist_ok=True)

    measurements = []
    for n in args.n_qubits:
        print(f"[orbit-gsim] n_qubits={n} weight={args.weight} gsim={gsim is not None}",
              file=sys.stderr, flush=True)
        measurements.append(
            measure_point(n, args.repeats, args.warmup, args.seed, args.weight, gsim)
        )

    payload = {
        "schema_version": 2,
        "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "label": args.label,
        "benchmark": "pauli_orbit_commutator_vs_gsim",
        "config": {
            "n_qubits_values": args.n_qubits,
            "weight": args.weight,
            "repeats": args.repeats,
            "warmup": args.warmup,
            "seed": args.seed,
            "gsim": gsim is not None,
        },
        "hardware": hardware.collect(),
        "measurements": measurements,
    }
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    png_path = output.with_suffix(".png")
    plot(payload, png_path)

    print(f"wrote {len(measurements)} measurements → {output}", file=sys.stderr)
    print(f"wrote plot → {png_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
