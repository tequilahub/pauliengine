"""Benchmark: PauliCycle commutator vs. the primitive QubitHamiltonian route.

A PauliCycle of ``n`` qubits stands for a base Pauli string and all ``n`` of its
cyclic rotations. There are two ways to obtain the commutator of two such cyclic
operators as a QubitHamiltonian:

* **cycle**  — stay in the Pauli-cycle basis: ``a.commutator(b)`` returns a
  PauliCycleSum without ever materialising the full operator. With bounded Pauli
  weight this costs O(n) (Corollary 3 of arXiv:2604.16701).

* **naive**  — expand first: build the full Hamiltonians ``H_a`` and ``H_b`` from
  all ``n`` rotations of each base, then ``H_a.commutator(H_b)``. This commutes
  every rotation of ``a`` against every rotation of ``b``.

Both must yield the identical operator (checked once per point). The sweep grows
the number of qubits ``n`` (which is also the number of rotations), so the plot
shows how the two approaches scale.

Usage (from the repo root):

    python -m tools.benchmarks.cycle_benchmark \
        --n-qubits 4 8 16 32 64 128 256 \
        --repeats 5 --warmup 1 \
        --output tools/benchmarks/results/cycle-latest.json

The PNG is written next to the JSON.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import random
import sys
from pathlib import Path

import pauliengine as pe
from pauliengine._core import PauliCycleComplex

# Support both invocation styles:
#   python -m tools.benchmarks.cycle_benchmark   (run as a package module)
#   python cycle_benchmark.py                    (run as a plain script)
if __package__:
    from . import hardware
    from .benchmark import _time_call
else:
    # Plain-script mode: put the repo root on sys.path and import via the
    # package, so intra-package imports inside benchmark.py resolve too.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from tools.benchmarks import hardware
    from tools.benchmarks.benchmark import _time_call

_PAULI_OPS = ("X", "Y", "Z")



def random_base(n_qubits: int, rng: random.Random, weight: int | None = None) -> str:
    """A random Pauli string of length ``n_qubits``.

    With ``weight`` set, exactly ``weight`` qubits carry a (random) non-identity
    operator and the rest are identities — this is the bounded-weight regime of
    Corollary 3 (e.g. the weight-1 / weight-2 generators of a TFIM ansatz),
    where the cycle commutator is O(n). ``weight=None`` (or weight >= n_qubits)
    gives a dense string, for which Corollary 3 offers no asymptotic saving.
    """
    if weight is None or weight >= n_qubits:
        return "".join(rng.choice(_PAULI_OPS) for _ in range(n_qubits))
    chars = ["I"] * n_qubits
    for p in rng.sample(range(n_qubits), weight):
        chars[p] = rng.choice(_PAULI_OPS)
    return "".join(chars)


def _rotations(base: str) -> list[str]:
    """All ``len(base)`` cyclic rotations: new[(i + k) mod n] = old[i]."""
    n = len(base)
    out = []
    for k in range(n):
        rot = [""] * n
        for i, ch in enumerate(base):
            rot[(i + k) % n] = ch
        out.append("".join(rot))
    return out


def _hamiltonian_from_base(base: str):
    """Primitive route: a QubitHamiltonian summing all rotations of ``base``."""
    terms = []
    for s in _rotations(base):
        ops = {i: ch for i, ch in enumerate(s) if ch != "I"}
        terms.append(pe.PauliString(complex(1.0), ops))
    return pe.QubitHamiltonian(terms)


def _make_naive_call(base_a: str, base_b: str):
    h_a = _hamiltonian_from_base(base_a)
    h_b = _hamiltonian_from_base(base_b)
    return lambda: h_a.commutator(h_b)


def measure_point(n_qubits: int, repeats: int, warmup: int, seed: int,
                  weight: int | None = None) -> dict:
    """Time the two routes at one qubit count and assert they agree.

    * ``naive``  — build both full Hamiltonians, then commute.
    * ``cycle``  — commutator in the Pauli-cycle basis, staying a PauliCycleSum
      (Eq. (40)). This is the fair like-for-like comparison against ``naive``:
      the same operator, in the symmetry-adapted basis. With bounded ``weight``
      this is O(n) (Corollary 3).
    """
    rng = random.Random(seed + n_qubits)
    base_a = random_base(n_qubits, rng, weight)
    base_b = random_base(n_qubits, rng, weight)

    a = PauliCycleComplex(n_qubits, base_a)
    b = PauliCycleComplex(n_qubits, base_b)

    naive_call = _make_naive_call(base_a, base_b)
    cycle_call = lambda: a.commutator(b)                          # stays compressed

    # Correctness check (not timed): materialising the cycle result must equal
    # the naive operator.
    if a.commutator(b).to_qubit_hamiltonian() != naive_call():
        raise AssertionError(f"cycle vs naive mismatch at n_qubits={n_qubits}")

    n_terms = len(naive_call())  # size of the resulting commutator operator

    n_min, n_mean, n_std = _time_call(naive_call, repeats=repeats, warmup=warmup)
    c_min, c_mean, c_std = _time_call(cycle_call, repeats=repeats, warmup=warmup)

    return {
        "n_qubits": n_qubits,
        "result_terms": n_terms,
        "repeats": repeats,
        "naive_time_min_s": n_min,
        "naive_time_mean_s": n_mean,
        "naive_time_stdev_s": n_std,
        "cycle_time_min_s": c_min,
        "cycle_time_mean_s": c_mean,
        "cycle_time_stdev_s": c_std,
        "speedup_mean": (n_mean / c_mean) if c_mean > 0 else float("nan"),
    }


def plot(payload: dict, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    ms = sorted(payload["measurements"], key=lambda m: m["n_qubits"])
    xs = [m["n_qubits"] for m in ms]
    weight = payload.get("config", {}).get("weight")

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: absolute runtime, log-log.
    for key, label, marker in (
        ("naive", "naive: full QubitHamiltonian commutator", "s"),
        ("cycle", f"cycle: PauliCycle", "o"),
    ):
        ys = [m[f"{key}_time_mean_s"] for m in ms]
        yerr = [m[f"{key}_time_stdev_s"] for m in ms]
        ax.errorbar(xs, ys, yerr=yerr, marker=marker, capsize=3, label=label)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    weight_label = "dense" if weight is None else f"weight {weight}"
    ax.set_xlabel("n_qubits  (= number of rotations)")
    ax.set_ylabel("commutator time [s] (mean ± stdev)")
    ax.set_title(f"Pauli-cycle vs. full QubitHamiltonian commutator ({weight_label})")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(fontsize=8)

    # Right: speedup factor.
    speedups = [m["speedup_mean"] for m in ms]
    ax2.plot(xs, speedups, marker="o", color="tab:green")
    ax2.axhline(1.0, color="gray", ls="--", alpha=0.6)
    ax2.set_xscale("log", base=2)
    ax2.set_xlabel("n_qubits")
    ax2.set_ylabel("speedup  (naive / cycle)")
    ax2.set_title("Speedup of the cycle route")
    ax2.grid(True, which="both", ls="--", alpha=0.4)

    hw_label = hardware.short_label(payload.get("hardware", {}))
    fig.text(0.5, 0.01, hw_label, ha="center", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="tools.benchmarks.cycle_benchmark",
        description="Benchmark PauliCycle commutator against the primitive QubitHamiltonian route.",
    )
    p.add_argument(
        "--n-qubits", type=int,
        default= 10,
        help="Qubit counts to sweep (also the number of rotations per cycle).",
    )
    p.add_argument(
        "--weight", type=int, default=3,
        help="Non-identity operators per base string (bounded-weight regime of "
             "Corollary 3). Use 0 for dense strings (weight = n_qubits).",
    )
    p.add_argument("--repeats", type=int, default=5, help="Measured runs per point.")
    p.add_argument("--warmup", type=int, default=1, help="Warm-up runs per point (untimed).")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--output", type=Path, default=None, help="Output JSON path.")
    p.add_argument("--label", type=str, default=None)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    output = args.output or (Path(__file__).parent / "results" / f"cycle-{stamp}.json")
    output.parent.mkdir(parents=True, exist_ok=True)

    weight = None if args.weight <= 0 else args.weight
    measurements = []
    for n in range(args.n_qubits):
        print(f"[cycle-bench] n_qubits={2**n} weight={weight}", file=sys.stderr, flush=True)
        measurements.append(
            measure_point(2**n, args.repeats, args.warmup, args.seed, weight)
        )

    payload = {
        "schema_version": 1,
        "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "label": args.label,
        "benchmark": "pauli_cycle_commutator",
        "config": {
            "n_qubits_values": args.n_qubits,
            "weight": weight,
            "repeats": args.repeats,
            "warmup": args.warmup,
            "seed": args.seed,
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
