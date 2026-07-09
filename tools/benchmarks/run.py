"""CLI entry point: run the benchmark sweeps and write a result JSON.

Usage (from the repo root):

    python -m tools.benchmarks.run \
        --n-terms 10 50 100 500 1000 \
        --n-qubits 4 8 16 32 64 \
        --n-terms-fixed 100 \
        --n-qubits-fixed 16 \
        --repeats 5 \
        --output tools/benchmarks/results/latest.json

Then plot with:

    python -m tools.benchmarks.plot tools/benchmarks/results/latest.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

from . import benchmark, hardware
from .generate import CoeffKind


def _default_output_path() -> Path:
    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    return Path(__file__).parent / "results" / f"bench-{stamp}.json"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="tools.benchmarks.run",
        description="Benchmark PauliEngine multiply and commutator on Complex + Symbolic Hamiltonians.",
    )
    p.add_argument(
        "--n-terms", type=int, nargs="+", default=[10, 50, 100, 500, 1000],
        help="Values swept for the n_terms axis (n_qubits held fixed).",
    )
    p.add_argument(
        "--n-qubits", type=int, nargs="+", default=[50, 100, 150, 200, 250],
        help="Values swept for the n_qubits axis (n_terms held fixed).",
    )
    p.add_argument(
        "--n-terms-fixed", type=int, default=200,
        help="Value of n_terms used during the n_qubits sweep.",
    )
    p.add_argument(
        "--n-qubits-fixed", type=int, default=100,
        help="Value of n_qubits used during the n_terms sweep.",
    )
    p.add_argument(
        "--ops", nargs="+", choices=["multiply", "commutator"],
        default=["multiply", "commutator"],
    )
    p.add_argument(
        "--coeff-kinds", nargs="+", choices=["complex", "symbolic"],
        default=["complex", "symbolic"],
    )
    p.add_argument("--repeats", type=int, default=5, help="Measured runs per point.")
    p.add_argument("--warmup", type=int, default=1, help="Warm-up runs per point (untimed).")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument(
        "--output", type=Path, default=None,
        help="Output JSON path. Default: results/bench-<timestamp>.json",
    )
    p.add_argument(
        "--label", type=str, default=None,
        help="Optional label stored in the result file, e.g. 'baseline' or 'PR-42'.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output = args.output or _default_output_path()
    output.parent.mkdir(parents=True, exist_ok=True)

    hw = hardware.collect()

    def progress(msg: str) -> None:
        print(msg, file=sys.stderr, flush=True)

    coeff_kinds: list[CoeffKind] = args.coeff_kinds
    results: list[benchmark.Measurement] = []

    results.extend(benchmark.scan_n_terms(
        ops=args.ops,
        coeff_kinds=coeff_kinds,
        n_terms_values=args.n_terms,
        n_qubits_fixed=args.n_qubits_fixed,
        repeats=args.repeats,
        warmup=args.warmup,
        seed=args.seed,
        progress=progress,
    ))
    results.extend(benchmark.scan_n_qubits(
        ops=args.ops,
        coeff_kinds=coeff_kinds,
        n_qubits_values=args.n_qubits,
        n_terms_fixed=args.n_terms_fixed,
        repeats=args.repeats,
        warmup=args.warmup,
        seed=args.seed,
        progress=progress,
    ))

    payload = {
        "schema_version": 1,
        "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "label": args.label,
        "config": {
            "n_terms_values": args.n_terms,
            "n_qubits_values": args.n_qubits,
            "n_terms_fixed": args.n_terms_fixed,
            "n_qubits_fixed": args.n_qubits_fixed,
            "ops": args.ops,
            "coeff_kinds": args.coeff_kinds,
            "repeats": args.repeats,
            "warmup": args.warmup,
            "seed": args.seed,
        },
        "hardware": hw,
        "measurements": [m.to_dict() for m in results],
    }
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {len(results)} measurements → {output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
