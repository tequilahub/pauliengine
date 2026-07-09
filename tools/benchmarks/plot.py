"""Plot benchmark results.

Generates one figure per (scaling_axis × operation). Each figure has one
curve per coefficient kind (complex / symbolic), on a log-log scale, with
error bars from the stdev across repeats. A footer stamps the hardware
label so plots from different machines aren't confused.

Usage:

    python -m tools.benchmarks.plot <result.json> [--out-dir tools/benchmarks/results]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

from . import hardware


def _group(measurements: list[dict]) -> dict[tuple[str, str], list[dict]]:
    """Group by (scaling_axis, op). Each group is one figure."""
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for m in measurements:
        groups[(m["scaling_axis"], m["op"])].append(m)
    return groups


def _plot_group(
    scaling_axis: str,
    op: str,
    entries: list[dict],
    hw_label: str,
    out_path: Path,
) -> None:
    by_kind: dict[str, list[dict]] = defaultdict(list)
    for e in entries:
        by_kind[e["coeff_kind"]].append(e)

    fig, ax = plt.subplots(figsize=(7, 5))
    for kind, points in sorted(by_kind.items()):
        points = sorted(points, key=lambda p: p[scaling_axis])
        xs = [p[scaling_axis] for p in points]
        ys = [p["time_mean_s"] for p in points]
        yerr = [p["time_stdev_s"] for p in points]
        ax.errorbar(xs, ys, yerr=yerr, marker="o", capsize=3, label=kind)

    ax.set_xlabel(scaling_axis)
    ax.set_ylabel("time [s] (mean ± stdev)")
    fixed_axis = "n_qubits" if scaling_axis == "n_terms" else "n_terms"
    fixed_value = entries[0][fixed_axis]
    ax.set_title(f"{op}  —  scaling with {scaling_axis}  ({fixed_axis}={fixed_value})")
    ax.grid(True, ls="--", alpha=0.4)
    ax.legend(title="coefficients")

    fig.text(0.5, 0.01, hw_label, ha="center", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tools.benchmarks.plot")
    parser.add_argument("result_json", type=Path, help="Path to a benchmark result JSON.")
    parser.add_argument(
        "--out-dir", type=Path, default=None,
        help="Directory for output PNGs. Default: alongside the result JSON.",
    )
    args = parser.parse_args(argv)

    if not args.result_json.exists():
        print(f"no such file: {args.result_json}", file=sys.stderr)
        return 1

    payload = json.loads(args.result_json.read_text(encoding="utf-8"))
    measurements = payload["measurements"]
    hw_label = hardware.short_label(payload.get("hardware", {}))

    out_dir = args.out_dir or args.result_json.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    stem = args.result_json.stem
    written: list[Path] = []
    for (scaling_axis, op), entries in _group(measurements).items():
        out_path = out_dir / f"{stem}__{op}__scaling-{scaling_axis}.png"
        _plot_group(scaling_axis, op, entries, hw_label, out_path)
        written.append(out_path)

    for p in written:
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
