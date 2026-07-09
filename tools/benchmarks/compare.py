"""Compare two or more benchmark result JSONs.

Produces:

* One figure per (scaling_axis, op, coeff_kind), with one curve per input file.
* A textual speedup table on stdout (relative to the first file = baseline).

If two runs come from different
CPUs, OSes, or pauliengine versions, a warning is printed. The plot footer
still shows every host so the difference is visible.

Usage:

    python -m tools.benchmarks.compare \
        results/baseline.json results/patched.json \
        [--out-dir results/comparison]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

from . import hardware


def _load(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["_path"] = path
    payload["_display"] = payload.get("label") or path.stem
    return payload


def _hardware_key(hw: dict) -> tuple:
    """Fingerprint used to decide if two runs are comparable."""
    cpu = (hw.get("cpu") or {}).get("model")
    os_name = (hw.get("os") or {}).get("system")
    pe_ver = (hw.get("packages") or {}).get("pauliengine")
    return (cpu, os_name, pe_ver)


def _warn_if_incomparable(runs: list[dict]) -> None:
    keys = {r["_display"]: _hardware_key(r.get("hardware", {})) for r in runs}
    unique = set(keys.values())
    if len(unique) > 1:
        print("WARNING: runs differ on (CPU, OS, pauliengine version):", file=sys.stderr)
        for name, key in keys.items():
            print(f"  {name}: cpu={key[0]!r} os={key[1]!r} pauliengine={key[2]!r}",
                  file=sys.stderr)
        print("Comparison may be misleading. Plot footer shows all hosts.",
              file=sys.stderr)


def _index(measurements: list[dict]) -> dict[tuple, dict]:
    """Index measurements by (scaling_axis, op, coeff_kind, n_terms, n_qubits)."""
    out: dict[tuple, dict] = {}
    for m in measurements:
        key = (m["scaling_axis"], m["op"], m["coeff_kind"], m["n_terms"], m["n_qubits"])
        out[key] = m
    return out


def _group_points(
    run: dict, scaling_axis: str, op: str, coeff_kind: str,
) -> list[dict]:
    matched = [
        m for m in run["measurements"]
        if m["scaling_axis"] == scaling_axis
        and m["op"] == op
        and m["coeff_kind"] == coeff_kind
    ]
    return sorted(matched, key=lambda p: p[scaling_axis])


def _all_combinations(runs: list[dict]) -> set[tuple[str, str, str]]:
    combos: set[tuple[str, str, str]] = set()
    for r in runs:
        for m in r["measurements"]:
            combos.add((m["scaling_axis"], m["op"], m["coeff_kind"]))
    return combos


def _plot_one(
    scaling_axis: str,
    op: str,
    coeff_kind: str,
    runs: list[dict],
    out_path: Path,
) -> bool:
    fig, ax = plt.subplots(figsize=(7, 5))
    plotted_any = False
    fixed_axis = "n_qubits" if scaling_axis == "n_terms" else "n_terms"
    fixed_values: set[int] = set()

    for run in runs:
        points = _group_points(run, scaling_axis, op, coeff_kind)
        if not points:
            continue
        plotted_any = True
        xs = [p[scaling_axis] for p in points]
        ys = [p["time_mean_s"] for p in points]
        yerr = [p["time_stdev_s"] for p in points]
        fixed_values.update(p[fixed_axis] for p in points)
        ax.errorbar(xs, ys, yerr=yerr, marker="o", capsize=3, label=run["_display"])

    if not plotted_any:
        plt.close(fig)
        return False

    ax.set_xlabel(scaling_axis)
    ax.set_ylabel("time [s] (mean ± stdev)")
    fixed_str = ",".join(str(v) for v in sorted(fixed_values))
    ax.set_title(f"{op} · {coeff_kind}  —  scaling with {scaling_axis}  ({fixed_axis}={fixed_str})")
    ax.grid(True, ls="--", alpha=0.4)
    ax.legend(title="run")

    footer = " || ".join(
        f"{r['_display']}: {hardware.short_label(r.get('hardware', {}))}" for r in runs
    )
    fig.text(0.5, 0.01, footer, ha="center", va="bottom", fontsize=6, color="gray", wrap=True)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return True


def _print_speedup_table(runs: list[dict]) -> None:
    """Print baseline-relative speedups for every measurement point.

    Speedup > 1 means the compared run is faster than the baseline at that point.
    """
    if len(runs) < 2:
        return

    baseline = runs[0]
    base_index = _index(baseline["measurements"])
    base_name = baseline["_display"]

    for other in runs[1:]:
        print()
        print(f"# speedup: {other['_display']}  vs  baseline={base_name}")
        print(f"{'scales_with':<12} {'op':<11} {'kind':<9} "
              f"{'n_terms':>7} {'n_qubits':>8} "
              f"{'base [s]':>10} {'this [s]':>10} {'speedup':>8}")
        for m in other["measurements"]:
            key = (m["scaling_axis"], m["op"], m["coeff_kind"], m["n_terms"], m["n_qubits"])
            base = base_index.get(key)
            if base is None:
                continue
            base_t = base["time_mean_s"]
            this_t = m["time_mean_s"]
            speedup = base_t / this_t if this_t > 0 else float("inf")
            marker = ""
            if speedup >= 1.05:
                marker = "  faster"
            elif speedup <= 0.95:
                marker = "  slower"
            print(f"{m['scaling_axis']:<12} {m['op']:<11} {m['coeff_kind']:<9} "
                  f"{m['n_terms']:>7} {m['n_qubits']:>8} "
                  f"{base_t:>10.4g} {this_t:>10.4g} {speedup:>7.2f}x{marker}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tools.benchmarks.compare")
    parser.add_argument("result_jsons", type=Path, nargs="+",
                        help="Two or more benchmark JSONs. First is the baseline.")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Directory for comparison PNGs. Default: alongside the first JSON.")
    parser.add_argument("--tag", type=str, default="compare",
                        help="Filename prefix for the generated PNGs.")
    args = parser.parse_args(argv)

    if len(args.result_jsons) < 2:
        print("need at least two JSON files to compare", file=sys.stderr)
        return 1

    for p in args.result_jsons:
        if not p.exists():
            print(f"no such file: {p}", file=sys.stderr)
            return 1

    runs = [_load(p) for p in args.result_jsons]
    _warn_if_incomparable(runs)

    out_dir = args.out_dir or args.result_jsons[0].parent
    out_dir.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for scaling_axis, op, coeff_kind in sorted(_all_combinations(runs)):
        out_path = out_dir / f"{args.tag}__{op}__{coeff_kind}__scaling-{scaling_axis}.png"
        if _plot_one(scaling_axis, op, coeff_kind, runs, out_path):
            written.append(out_path)

    for p in written:
        print(p)

    _print_speedup_table(runs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
