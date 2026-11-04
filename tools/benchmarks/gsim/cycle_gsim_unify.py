"""Cost of unifying Pauli-cycle results: PauliCycle commutator / product with and
without ``PauliCycleSum.unify()``, compared with g-sim and the naive
QubitHamiltonian route, up to n = 2^12.

Routes per operation:

* **naive**        — expand both cycles into QubitHamiltonians and commute / multiply
* **cycle**        — ``PauliCycle.commutator`` / ``*`` (stored bases as produced)
* **cycle+unify**  — the same, followed by ``PauliCycleSum.unify()`` so every summand
  carries the unified representative of its orbit
* **gsim**         — g-sim's ``cycle_commutator`` (commutator only; g-sim has no cycle
  product). g-sim always returns canonical cycles, so *cycle+unify* is the
  like-for-like comparison.

Expensive routes are capped at a smaller n (their curves end there): the naive
route materialises up to n^2 strings, and dense g-sim grows ~n^3. Repeats adapt to
the per-call cost. Raw samples are stored (``*_times_s``); power laws are fitted on
n >= 32. Each figure has the runtimes on top and the overhead of unify below.

Run (from the repo root, g-sim optional):

    python -m tools.benchmarks.gsim.cycle_gsim_unify
"""

from __future__ import annotations

import datetime as dt
import json
import random
import sys
import time
from pathlib import Path

import numpy as np

if __package__:
    from ..utils import hardware
    from ..utils.benchmark import _time_call, route_stats
    from ..utils.gsim_compat import load_cycles
    from ..native.cycle_benchmark import random_base, _make_naive_call, _hamiltonian_from_base
    from .cycle_gsim_benchmark import _base_to_vec, _gsim_result_to_qh
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from tools.benchmarks.utils import hardware
    from tools.benchmarks.utils.benchmark import _time_call, route_stats
    from tools.benchmarks.utils.gsim_compat import load_cycles
    from tools.benchmarks.native.cycle_benchmark import random_base, _make_naive_call, _hamiltonian_from_base
    from tools.benchmarks.gsim.cycle_gsim_benchmark import _base_to_vec, _gsim_result_to_qh

from pauliengine import PauliCycle

RESULTS = Path(__file__).resolve().parents[1] / "results"
EXPONENTS = list(range(2, 13))          # n = 4 .. 4096
FIT_FROM = 32
BUDGET_S = 4.0                          # target measuring time per route and point
GSIM_CHECK_MAX_N = 64

# (operation, weight label, weight, caps for the expensive routes)
RUNS = [
    ("commutator", "weight3", 3, {"naive": 4096, "gsim": 4096}),
    ("commutator", "dense", None, {"naive": 1024, "gsim": 1024}),
    ("multiply", "weight3", 3, {"naive": 512, "gsim": 0}),
    ("multiply", "dense", None, {"naive": 512, "gsim": 0}),
]
ROUTES = [
    ("naive", "QubitHamiltonian (naive)", "tab:blue", "s"),
    ("cycle", "PauliCycle", "tab:orange", "o"),
    ("cycle_unify", "PauliCycle + unify()", "tab:red", "D"),
    ("gsim", "g-sim cycle_commutator", "tab:green", "^"),
]


def adaptive_times(fn) -> list[float]:
    t0 = time.perf_counter()
    fn()                                # also serves as warm-up
    one = time.perf_counter() - t0
    repeats = int(min(15, max(3, BUDGET_S / max(one, 1e-9))))
    return _time_call(fn, repeats=repeats, warmup=2 if one < 0.05 else 0)


def sweep(op: str, weight, caps: dict, gsim) -> dict:
    measurements = []
    for e in EXPONENTS:
        n = 2 ** e
        rng = random.Random(n)
        base_a, base_b = random_base(n, rng, weight), random_base(n, rng, weight)
        a, b = PauliCycle(n, base_a), PauliCycle(n, base_b)
        cyc = (lambda: a.commutator(b)) if op == "commutator" else (lambda: a * b)
        rec = {"n_qubits": n}

        if n <= caps["naive"]:
            if op == "commutator":
                naive = _make_naive_call(base_a, base_b)
            else:
                h_a, h_b = _hamiltonian_from_base(base_a), _hamiltonian_from_base(base_b)
                naive = lambda: h_a * h_b
            ref = naive()
            unified = cyc()
            unified.unify()
            if cyc().to_qubit_hamiltonian() != ref or unified.to_qubit_hamiltonian() != ref:
                raise AssertionError(f"{op}: cycle vs naive mismatch at n={n}")
            rec["naive_times_s"] = adaptive_times(naive)

        rec["cycle_times_s"] = adaptive_times(cyc)
        rec["cycle_unify_times_s"] = adaptive_times(lambda: cyc().unify())

        if op == "commutator" and gsim is not None and n <= caps["gsim"]:
            cc, canon = gsim
            vp, vq, cache = _base_to_vec(base_a, n, canon), _base_to_vec(base_b, n, canon), {}
            g = lambda: cc(vp, vq, n, cache)
            if n <= min(GSIM_CHECK_MAX_N, caps["naive"]) and _gsim_result_to_qh(g(), n) != ref:
                raise AssertionError(f"gsim vs naive mismatch at n={n}")
            rec["gsim_times_s"] = adaptive_times(g)

        routes = [k[:-8] for k in rec if k.endswith("_times_s")]
        print(f"[{op} {'dense' if weight is None else f'w{weight}'}] n={n:5}  "
              + "  ".join(f"{r}={route_stats(rec, r)[0] * 1e3:10.4f}ms" for r in routes), flush=True)
        measurements.append(rec)
    return {
        "schema_version": 2,
        "timestamp": dt.datetime.now(dt.timezone.utc).isoformat(),
        "benchmark": f"pauli_cycle_{op}_unify",
        "config": {"op": op, "weight": weight, "exponents": EXPONENTS, "caps": caps,
                   "budget_s": BUDGET_S, "fit_from": FIT_FROM, "gsim": gsim is not None},
        "hardware": hardware.collect(),
        "measurements": measurements,
    }


def power_fit(n: np.ndarray, y: np.ndarray):
    mask = n >= FIT_FROM
    if mask.sum() < 2:
        return None
    x, yy = np.log(n[mask]), np.log(y[mask])
    alpha, b = np.polyfit(x, yy, 1)
    r2 = 1 - np.sum((yy - (alpha * x + b)) ** 2) / np.sum((yy - yy.mean()) ** 2)
    return alpha, float(np.exp(b)), r2, n[mask]


def plot(payload: dict, out: Path) -> dict:
    import matplotlib.pyplot as plt

    cfg = payload["config"]
    ms = payload["measurements"]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(9, 8.5), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1.2]})
    fits = {}
    for key, label, color, marker in ROUTES:
        pts = [m for m in ms if f"{key}_times_s" in m]
        if not pts:
            continue
        n = np.array([m["n_qubits"] for m in pts], float)
        y, yerr = (np.array(v, float) for v in zip(*(route_stats(m, key) for m in pts)))
        f = power_fit(n, y)
        ax.errorbar(n, y, yerr=yerr, marker=marker, capsize=3, color=color,
                    label=label + (f"   α={f[0]:.2f}" if f else ""))
        if f:
            fits[key] = {"alpha": f[0], "C": f[1], "r2": f[2], "n_max": int(n.max())}
            xf = np.array([f[3].min(), f[3].max()])
            ax.plot(xf, f[1] * xf ** f[0], ls="--", color=color, alpha=0.6)

    wl = "dense" if cfg["weight"] is None else f"weight {cfg['weight']}"
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_ylabel(f"{cfg['op']} time [s] (geo. mean ± geo. spread)")
    ax.set_title(f"Pauli-cycle {cfg['op']} with / without unify() — {wl}  (fit on n ≥ {FIT_FROM})")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(fontsize=8, loc="upper left")
    caps = {k: v for k, v in cfg["caps"].items() if v}
    note = ", ".join(f"{k} to n={v}" for k, v in caps.items() if v < 2 ** EXPONENTS[-1])
    if cfg["op"] == "multiply":
        note = (note + "; " if note else "") + "g-sim has no cycle product"
    if note:
        ax.text(0.98, 0.03, "measured: " + note, transform=ax.transAxes, ha="right", va="bottom",
                fontsize=7, color="dimgray")

    # overhead of unify relative to the plain cycle route
    n = np.array([m["n_qubits"] for m in ms], float)
    base = np.array([route_stats(m, "cycle")[0] for m in ms])
    uni = np.array([route_stats(m, "cycle_unify")[0] for m in ms])
    over = (uni / base - 1.0) * 100.0
    ax2.bar(n, over, width=n * 0.35, color="tab:red", alpha=0.7)
    ax2.axhline(0, color="gray", lw=0.8)
    ax2.set_xscale("log", base=2)
    ax2.set_xlabel("n_qubits  (= number of rotations)")
    ax2.set_ylabel("unify overhead [%]")
    ax2.grid(True, axis="y", ls="--", alpha=0.4)

    fig.text(0.5, 0.005, hardware.short_label(payload.get("hardware", {})),
             ha="center", va="bottom", fontsize=7, color="gray")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out, dpi=150)
    plt.close(fig)
    fits["unify_overhead_pct"] = {int(k): float(v) for k, v in zip(n, over)}
    return fits


def main() -> int:
    gsim = load_cycles()
    if gsim is None:
        print("[cycle-unify] g-sim not available — skipping the g-sim route.", file=sys.stderr)
    RESULTS.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d")
    summary = {}
    for op, wname, weight, caps in RUNS:
        payload = sweep(op, weight, caps, gsim)
        stem = RESULTS / f"cycle-unify-{op}-{wname}-{stamp}"
        stem.with_suffix(".json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        summary[f"{op} {wname}"] = plot(payload, stem.with_suffix(".png"))
        print(f"  -> {stem.with_suffix('.png')}", flush=True)

    print(f"\n=== Power-law exponents (fit on n >= {FIT_FROM}) and unify overhead ===")
    for name, fits in summary.items():
        print(name)
        for k, f in fits.items():
            if k == "unify_overhead_pct":
                print("   unify overhead [%]: " + ", ".join(f"n={n}: {v:+.1f}" for n, v in f.items()))
            else:
                print(f"   {k:12} alpha = {f['alpha']:5.2f}  R2 = {f['r2']:.3f}  (n <= {f['n_max']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
