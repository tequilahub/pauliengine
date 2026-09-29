"""Benchmark: Pauli-cycle commutator — three routes compared.

Same setup as ``cycle_benchmark`` but with g-sim added as a third route:

* **naive**  — expand each cycle into an explicit QubitHamiltonian (all n
  rotations of the base) and commute the two full operators.
* **cycle**  — PauliEngine's Pauli-cycle commutator (``PauliCycle.commutator``),
  in C++.
* **gsim**   — g-sim's ``cycle_commutator`` (Python, from
  https://github.com/adelina-b/g-sim, src/gsim/gcycles.py).

Note: g-sim's cycle routine is a pure-Python loop that itself calls a PauliEngine
``PauliString.commutator`` per active shift (``import PauliEngine as pe`` inside
gcycles), so this measures the native C++ cycle basis against a Python-level
driver over single-string commutators.

All three describe the same operator; correctness is cross-checked once per point
(orbit route vs the expanded operator always; g-sim vs the expanded operator for
small n only, to keep the check cheap).

g-sim is optional and used automatically if the ``gsim`` package is installed.
``gsim_compat`` registers a ``PauliEngine`` shim so g-sim's cycle module (which
imports ``PauliEngine`` and constructs strings as ``PauliString(dict, coeff)``)
runs against the locally installed ``pauliengine``. If unavailable, the g-sim
route is skipped.

Usage (from the repo root):

    python -m tools.benchmarks.cycle_gsim_benchmark \
        --n-qubits 4 8 16 32 64 128 --weight 3 \
        --repeats 5 --warmup 1 \
        --output tools/benchmarks/results/cycle-gsim-latest.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import random
import sys
from pathlib import Path

import pauliengine as pe
from pauliengine import PauliCycle

if __package__:
    from ..utils import hardware
    from ..utils.benchmark import _time_call, has_route, route_stats
    from ..utils.gsim_compat import load_cycles
    from ..native.cycle_benchmark import random_base, _rotations, _make_naive_call
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from tools.benchmarks.utils import hardware
    from tools.benchmarks.utils.benchmark import _time_call, has_route, route_stats
    from tools.benchmarks.utils.gsim_compat import load_cycles
    from tools.benchmarks.native.cycle_benchmark import random_base, _rotations, _make_naive_call

_P2I = {"X": 1, "Y": 2, "Z": 3}
_I2P = {1: "X", 2: "Y", 3: "Z"}

# Above this qubit count the (untimed) g-sim correctness check — which expands the
# result back to a full operator — is skipped to keep it cheap.
_GSIM_CHECK_MAX_N = 64


def _base_to_vec(base: str, n: int, get_canonical_cycle) -> dict:
    """A single cycle (base string) as g-sim's sparse vector {canonical_sig: 1.0}."""
    sig = tuple(sorted((i, _P2I[ch]) for i, ch in enumerate(base) if ch != "I"))
    return {get_canonical_cycle(sig, n): 1.0}


def _gsim_result_to_qh(result: dict, n: int):
    """Expand g-sim's {canonical_sig: coeff} result into a QubitHamiltonian."""
    terms = []
    for sig, coeff in result.items():
        base = ["I"] * n
        for idx, code in sig:
            base[idx] = _I2P[code]
        for rot in _rotations("".join(base)):
            ops = {i: ch for i, ch in enumerate(rot) if ch != "I"}
            terms.append(pe.PauliString(complex(coeff), ops))
    return pe.QubitHamiltonian(terms)


def measure_point(n_qubits: int, repeats: int, warmup: int, seed: int,
                  weight: int | None, gsim=None) -> dict:
    """Time the routes at one qubit count and assert they agree."""
    rng = random.Random(seed + n_qubits)
    base_a = random_base(n_qubits, rng, weight)
    base_b = random_base(n_qubits, rng, weight)

    a = PauliCycle(n_qubits, base_a)
    b = PauliCycle(n_qubits, base_b)

    naive_call = _make_naive_call(base_a, base_b)
    cycle_call = lambda: a.commutator(b)

    if a.commutator(b).to_qubit_hamiltonian() != naive_call():
        raise AssertionError(f"cycle vs naive mismatch at n_qubits={n_qubits}")

    record = {
        "n_qubits": n_qubits,
        "result_terms": len(naive_call()),
        "repeats": repeats,
    }

    for name, call in (("naive", naive_call), ("cycle", cycle_call)):
        record[f"{name}_times_s"] = _time_call(call, repeats=repeats, warmup=warmup)

    if gsim is not None:
        cycle_commutator, get_canonical_cycle = gsim
        vec_p = _base_to_vec(base_a, n_qubits, get_canonical_cycle)
        vec_q = _base_to_vec(base_b, n_qubits, get_canonical_cycle)
        ps_cache: dict = {}  # persistent across warmup+timed runs (g-sim's design)
        gsim_call = lambda: cycle_commutator(vec_p, vec_q, n_qubits, ps_cache)

        # Correctness (small n only): expand g-sim's result and compare operators.
        if n_qubits <= _GSIM_CHECK_MAX_N:
            if _gsim_result_to_qh(gsim_call(), n_qubits) != naive_call():
                raise AssertionError(f"gsim vs naive mismatch at n_qubits={n_qubits}")

        record["gsim_times_s"] = _time_call(gsim_call, repeats=repeats, warmup=warmup)

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
        ("cycle", "cycle: PauliEngine PauliCycle", "o", "tab:orange"),
    ]
    if has_gsim:
        routes.append(("gsim", "gsim: g-sim cycle_commutator (python)", "^", "tab:green"))

    for key, label, marker, color in routes:
        ys, yerr = zip(*(route_stats(m, key) for m in ms))
        ax.errorbar(xs, ys, yerr=yerr, marker=marker, capsize=3, color=color, label=label)

    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("n_qubits  (= number of rotations)")
    ax.set_ylabel("commutator time [s] (geo. mean ± geo. spread)")
    weight_label = "dense" if weight is None else f"weight {weight}"
    ax.set_title(f"Pauli-cycle commutator: PauliEngine vs g-sim vs naive ({weight_label})")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(fontsize=8)

    hw_label = hardware.short_label(payload.get("hardware", {}))
    fig.text(0.5, 0.01, hw_label, ha="center", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="tools.benchmarks.cycle_gsim_benchmark",
        description="Benchmark the Pauli-cycle commutator: PauliEngine vs g-sim vs naive.",
    )
    p.add_argument("--n-qubits", type=int, nargs="+", default=[4, 8, 16, 32, 64, 128])
    p.add_argument("--weight", type=int, default=3,
                   help="Non-identity operators per base string (0 => dense).")
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--warmup", type=int, default=1)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--output", type=Path, default=None)
    p.add_argument("--label", type=str, default=None)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    gsim = load_cycles()
    if gsim is None:
        print("[cycle-gsim] g-sim not available — measuring naive + cycle only.",
              file=sys.stderr)

    weight = None if args.weight <= 0 else args.weight

    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    output = args.output or (Path(__file__).resolve().parents[1] / "results" / f"cycle-gsim-{stamp}.json")
    output.parent.mkdir(parents=True, exist_ok=True)

    measurements = []
    for n in args.n_qubits:
        print(f"[cycle-gsim] n_qubits={n} weight={weight} gsim={gsim is not None}",
              file=sys.stderr, flush=True)
        measurements.append(
            measure_point(n, args.repeats, args.warmup, args.seed, weight, gsim)
        )

    payload = {
        "schema_version": 2,
        "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "label": args.label,
        "benchmark": "pauli_cycle_commutator_vs_gsim",
        "config": {
            "n_qubits_values": args.n_qubits,
            "weight": weight,
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
