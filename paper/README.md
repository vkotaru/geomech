# Scripts of the geomech paper

These scripts compute the numbers that the paper about geomech reports. Run them from the
repository root with geomech installed (`python -m pip install -e ".[dev]"`; SciPy is needed
for the control examples).

| Script | What it does |
|---|---|
| `sim/check_appendix.py` | Compares every equation displayed in the catalog appendix with the output of `compute_eom`, at default and random parameters; fails above 1e-12 relative error. |
| `sim/check_references.py` | Compares the quadrotor with Lee et al. (2010) and the flexible cable (n = 1, 3, 5) with Kotaru et al. (2018). |
| `sim/check_symmetry.py` | Checks that the standard-form mass matrix of every library system is symmetric along a short simulation. |
| `sim/catalog_numbers.py` | Derivation time and numerical check of every catalog entry (the verification table). Timings vary between runs and machines. |
| `sim/scaling_numbers.py` | Derivation time and output size against the number of bodies. |
| `sim/quadrotor_load_lqr.py` | The LQR listing of the quadrotor with a suspended load. |
| `sim/control/export.py` | Turns the stored results of `examples/control/run.py` into the control table, macros and plots. |

```bash
python paper/sim/check_appendix.py
python examples/control/run.py --all
python paper/sim/control/export.py examples/control/results
```

The scripts write their tables and data files next to themselves; those outputs are not tracked.
