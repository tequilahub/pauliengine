"""Benchmark: PauliOrbit commutator vs. the primitive QubitHamiltonian route.

A Pauli orbit ``B_{p,q,r}`` on ``n`` qubits stands for the (permutation-invariant)
sum of ALL Pauli strings carrying X, Y, Z on exactly p, q, r sites. There are two
ways to obtain the commutator of two such orbits:

* **orbit**  — stay in the Pauli-orbit basis: ``a.commutator(b)`` returns a
  PauliOrbitSum via the contingency-table formula, without ever materialising the
  exponentially many Pauli strings (Theorem 1 / Prop. 7 of arXiv:2604.16701). For
  bounded weight this is essentially constant per output coefficient (Cor. 5).

* **naive**  — expand first: build the full Hamiltonians ``H_a`` and ``H_b`` as the
  ``N_terms`` explicit Pauli strings of each orbit, then ``H_a.commutator(H_b)``.

Both must yield the identical operator (checked once per point, up to a tolerance
because the orbit coefficients are rationals). The sweep grows ``n``. Because a
bounded-weight orbit already expands to ``~n^weight`` Pauli strings — and the naive
commutator to ``~n^(2*weight)`` pair products — the qubit counts are kept much
SMALLER than in the Pauli-cycle benchmark.

Usage (from the repo root):

    python -m tools.benchmarks.orbit_benchmark \
        --n-qubits 4 6 8 12 16 24 32 \
        --weight 2 --repeats 5 --warmup 1 \
        --output tools/benchmarks/results/orbit-latest.json

The PNG is written next to the JSON.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import random
import sys
from itertools import combinations
from pathlib import Path

import pauliengine as pe
from pauliengine import PauliOrbit

# Support both invocation styles:
#   python -m tools.benchmarks.orbit_benchmark   (run as a package module)
#   python orbit_benchmark.py                    (run as a plain script)
if __package__:
    from ..utils import hardware
    from ..utils.benchmark import _time_call, route_stats
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from tools.benchmarks.utils import hardware
    from tools.benchmarks.utils.benchmark import _time_call, route_stats


def n_terms(n: int, p: int, q: int, r: int) -> int:
    """Number of distinct Pauli strings in orbit B_{p,q,r} (Eq. 52)."""
    return math.comb(n, p) * math.comb(n - p, q) * math.comb(n - p - q, r)


def random_label(n: int, rng: random.Random, weight: int) -> tuple[int, int, int]:
    """A random orbit label (p, q, r) with p + q + r == min(weight, n).

    Fixing the total weight keeps the naive expansion at ~n^weight strings, i.e.
    the bounded-weight regime where the orbit route is expected to win.
    """
    w = min(weight, n)
    cuts = sorted(rng.randint(0, w) for _ in range(2))
    return cuts[0], cuts[1] - cuts[0], w - cuts[1]


def orbit_to_hamiltonian(n: int, p: int, q: int, r: int):
    """Primitive route: B_{p,q,r} as (i / N_terms) * sum over all placements.

    Placements are enumerated directly via ``combinations`` (choosing the X, Y,
    then Z sites), which produces exactly ``N_terms`` strings — never the ``n!``
    that a naive ``permutations`` would.
    """
    N = n_terms(n, p, q, r)
    coeff = 1j / N
    terms = []
    all_sites = range(n)
    for x_pos in combinations(all_sites, p):
        x_set = set(x_pos)
        rem1 = [i for i in all_sites if i not in x_set]
        for y_pos in combinations(rem1, q):
            y_set = set(y_pos)
            rem2 = [i for i in rem1 if i not in y_set]
            for z_pos in combinations(rem2, r):
                ops = {i: "X" for i in x_pos}
                ops.update({i: "Y" for i in y_pos})
                ops.update({i: "Z" for i in z_pos})
                terms.append(pe.PauliString(coeff, ops))
    return pe.QubitHamiltonian(terms)


def _orbitsum_to_hamiltonian(s):
    """Expand a PauliOrbitSum back into a QubitHamiltonian (for the check)."""
    acc = pe.QubitHamiltonian.zero()
    for orbit, c in zip(s.data, s.coeffs):
        acc = acc + orbit_to_hamiltonian(s.n_qubits, orbit.p, orbit.q, orbit.r) * complex(c)
    return acc


def _qh_close(a, b, tol: float = 1e-9) -> bool:
    """Coefficient-wise comparison of two QubitHamiltonians within a tolerance."""

    def as_dict(qh):
        return {tuple(sorted(ops.items())): complex(c) for c, ops in qh.to_dictionary()}

    da, db = as_dict(a), as_dict(b)
    return all(abs(da.get(k, 0j) - db.get(k, 0j)) < tol for k in set(da) | set(db))


def measure_point(n_qubits: int, repeats: int, warmup: int, seed: int,
                  weight: int = 2) -> dict:
    """Time the two routes at one qubit count and assert they agree.

    * ``naive``  — build both full Hamiltonians, then commute.
    * ``orbit``  — commutator in the Pauli-orbit basis, staying a PauliOrbitSum.
    """
    rng = random.Random(seed + n_qubits)
    label_a = random_label(n_qubits, rng, weight)
    label_b = random_label(n_qubits, rng, weight)

    a = PauliOrbit(n_qubits, *label_a)
    b = PauliOrbit(n_qubits, *label_b)

    h_a = orbit_to_hamiltonian(n_qubits, *label_a)
    h_b = orbit_to_hamiltonian(n_qubits, *label_b)
    naive_call = lambda: h_a.commutator(h_b)
    orbit_call = lambda: a.commutator(b)                          # stays compressed

    # Correctness check (not timed): materialising the orbit result must equal
    # the naive operator (up to a tolerance; orbit coefficients are rationals).
    if not _qh_close(_orbitsum_to_hamiltonian(a.commutator(b)), naive_call()):
        raise AssertionError(f"orbit vs naive mismatch at n_qubits={n_qubits}")

    result_terms = len(naive_call())     # size of the resulting commutator operator

    return {
        "n_qubits": n_qubits,
        "label_a": list(label_a),
        "label_b": list(label_b),
        "orbit_strings": n_terms(n_qubits, *label_a),
        "result_terms": result_terms,
        "repeats": repeats,
        "naive_times_s": _time_call(naive_call, repeats=repeats, warmup=warmup),
        "orbit_times_s": _time_call(orbit_call, repeats=repeats, warmup=warmup),
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
        ("orbit", "orbit: PauliOrbit", "o"),
    ):
        ys, yerr = zip(*(route_stats(m, key) for m in ms))
        ax.errorbar(xs, ys, yerr=yerr, marker=marker, capsize=3, label=label)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("n_qubits")
    ax.set_ylabel("commutator time [s] (geo. mean ± geo. spread)")
    ax.set_title(f"Pauli-orbit vs. full QubitHamiltonian commutator (weight {weight})")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(fontsize=8)

    # Right: speedup factor.
    speedups = [route_stats(m, "naive")[0] / route_stats(m, "orbit")[0] for m in ms]
    ax2.plot(xs, speedups, marker="o", color="tab:green")
    ax2.axhline(1.0, color="gray", ls="--", alpha=0.6)
    ax2.set_xscale("log", base=2)
    ax2.set_xlabel("n_qubits")
    ax2.set_ylabel("speedup  (naive / orbit)")
    ax2.set_title("Speedup of the orbit route")
    ax2.grid(True, which="both", ls="--", alpha=0.4)

    hw_label = hardware.short_label(payload.get("hardware", {}))
    fig.text(0.5, 0.01, hw_label, ha="center", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="tools.benchmarks.orbit_benchmark",
        description="Benchmark PauliOrbit commutator against the primitive QubitHamiltonian route.",
    )
    p.add_argument(
        "--n-qubits", type=int, nargs="+",
        default=[4, 6, 8, 12, 16, 24, 32],
        help="Qubit counts to sweep. Kept small: orbits expand to ~n^weight strings.",
    )
    p.add_argument(
        "--weight", type=int, default=2,
        help="Total non-identity weight p+q+r of each orbit label (bounded-weight regime).",
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
    output = args.output or (Path(__file__).resolve().parents[1] / "results" / f"orbit-{stamp}.json")
    output.parent.mkdir(parents=True, exist_ok=True)

    measurements = []
    for n in args.n_qubits:
        print(f"[orbit-bench] n_qubits={n} weight={args.weight}", file=sys.stderr, flush=True)
        measurements.append(
            measure_point(n, args.repeats, args.warmup, args.seed, args.weight)
        )

    payload = {
        "schema_version": 2,
        "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "label": args.label,
        "benchmark": "pauli_orbit_commutator",
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
    plot(payload, png_path)

    print(f"wrote {len(measurements)} measurements → {output}", file=sys.stderr)
    print(f"wrote plot → {png_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
