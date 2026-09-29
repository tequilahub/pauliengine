# PauliEngine Benchmarks

Micro-benchmarks for `QubitHamiltonian * QubitHamiltonian` (multiply) and
`QubitHamiltonian.commutator(...)`, run against both the complex and
symbolic coefficient backends.

Each run produces one self-contained JSON: measurements plus the hardware
and package fingerprint of the machine that produced them. Never compare
plots across machines by eye — check the hardware footer first.

## Layout

```
tools/benchmarks/
├── utils/          # shared infrastructure (import-only, no CLIs of their own except plot/compare)
│   ├── hardware.py     # host + version fingerprint
│   ├── generate.py     # reproducible random Hamiltonians
│   ├── benchmark.py    # timing (geometric mean/spread) + scan primitives
│   ├── gsim_compat.py  # loads g-sim + PauliEngine import shim
│   ├── plot.py         # CLI, writes results/*.png next to the JSON
│   └── compare.py      # CLI, overlays several runs
├── native/         # PauliEngine benchmarks
│   ├── run.py             # QubitHamiltonian multiply/commutator scan (CLI)
│   ├── cycle_benchmark.py / cycle_fit_plot.py     # PauliCycle commutator
│   └── orbit_benchmark.py / orbit_fit_plot.py     # PauliOrbit commutator
├── gsim/           # PauliEngine-vs-g-sim comparison benchmarks
│   ├── cycle_gsim_benchmark.py / cycle_gsim_fit_plot.py
│   └── orbit_gsim_benchmark.py / orbit_gsim_fit_plot.py
└── results/        # outputs land here (shared across all benchmarks)
```

## Requirements

```
pip install matplotlib
```

PauliEngine itself must be importable — install it (`pip install -e .`
from the repo root) before running.

## Run a scan

From the repo root:

```bash
python -m tools.benchmarks.native.run \
    --n-terms 10 50 100 500 1000 \
    --n-qubits 4 8 16 32 \
    --n-terms-fixed 100 \
    --n-qubits-fixed 8 \
    --repeats 5
```

The two scan axes are independent:

* `--n-terms` scans the number of Pauli strings, keeping length fixed at `--n-qubits-fixed`.
* `--n-qubits` scans the length of each Pauli string, keeping count fixed at `--n-terms-fixed`.

Symbolic scans get expensive fast — for a quick check pass
`--coeff-kinds complex` and small scan values.

Add `--label baseline` (or similar) to tag the run in the JSON.

## Plot

```bash
python -m tools.benchmarks.utils.plot tools/benchmarks/results/bench-YYYYMMDD-HHMMSS.json
```

Writes one PNG per `(scaling_axis × op)` combination next to the JSON. Both
coefficient kinds appear as separate curves in each figure, linear axes,
with error bars from the stdev across repeats. The CPU model, OS, and
pauliengine version are stamped at the bottom of every plot.

## Compare against older runs

```bash
python -m tools.benchmarks.utils.compare \
    tools/benchmarks/results/baseline.json \
    tools/benchmarks/results/patched.json
```

The first file is the baseline. For every `(scaling_axis, op, coeff_kind)`
combination that appears in either file, one PNG is written with all runs
overlaid (linear axes). A speedup table is printed to stdout (values > 1× mean the later
run is faster than the baseline at that point).

If two runs come from different CPUs, OSes, or `pauliengine` versions, a
warning is printed and every host appears in the plot footer — the tool
never silently glues incomparable data together. Compare within the same
machine whenever possible.

Points that exist in one file but not another are simply skipped on the
speedup table; the plot still shows whichever runs have data.

## Notes

* Timing uses `time.perf_counter`, with GC disabled during measurement and
  a warm-up run to prime caches. The JSON stores **every raw sample**
  (`times_s`, or `<route>_times_s` for the cycle/orbit/g-sim benchmarks)
  instead of precomputed statistics; plots derive the geometric mean and
  spread from them (`utils.benchmark.route_stats`). Result files from before
  this change (`schema_version` 1, `*_time_mean_s`) can still be plotted.
* Seeds are fixed by default so two runs on the same machine hit the exact
  same Hamiltonians.
* The result schema is versioned (`schema_version` in the JSON) so `plot.py`
  can evolve without silently misreading old runs.

## Cycle / orbit / g-sim benchmarks

Beyond the QubitHamiltonian scan above, the symmetry-adapted primitives have
their own benchmarks (each writes JSON + PNG into `results/`):

```bash
# PauliEngine only: commutator in the symmetry-adapted basis vs. the naive
# expanded-QubitHamiltonian route, with a power-law fit.
python -m tools.benchmarks.native.cycle_fit_plot
python -m tools.benchmarks.native.orbit_fit_plot

# PauliEngine vs g-sim (needs the `gsim` package installed; loaded automatically
# via utils/gsim_compat.py, no path to configure).
python -m tools.benchmarks.gsim.cycle_gsim_fit_plot
python -m tools.benchmarks.gsim.orbit_gsim_fit_plot
```

Timing aggregates each point by the geometric mean of the samples with a
geometric-standard-deviation spread (`utils/benchmark.py`), which is more robust
for the positive, right-skewed runtimes.
