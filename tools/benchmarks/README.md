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
├── hardware.py     # host + version fingerprint
├── generate.py     # reproducible random Hamiltonians
├── benchmark.py    # timing + scan primitives
├── run.py          # CLI, writes results/*.json
├── plot.py         # CLI, writes results/*.png next to the JSON
└── results/        # outputs land here (gitignored by convention)
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
python -m tools.benchmarks.run \
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
python -m tools.benchmarks.plot tools/benchmarks/results/bench-YYYYMMDD-HHMMSS.json
```

Writes one PNG per `(scaling_axis × op)` combination next to the JSON. Both
coefficient kinds appear as separate curves in each figure, linear axes,
with error bars from the stdev across repeats. The CPU model, OS, and
pauliengine version are stamped at the bottom of every plot.

## Compare against older runs

```bash
python -m tools.benchmarks.compare \
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
  a warm-up run to prime caches. Report the *minimum* if you're chasing
  best-case throughput, the *mean* for typical behavior — both are stored.
* Seeds are fixed by default so two runs on the same machine hit the exact
  same Hamiltonians.
* The result schema is versioned (`schema_version` in the JSON) so `plot.py`
  can evolve without silently misreading old runs.
