# Contributing to geomech

Thanks for your interest in improving geomech. Bug reports, new verified systems,
documentation fixes and engine improvements are all welcome.

## Development setup

```bash
git clone https://github.com/vkotaru/geomech.git
cd geomech
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

## Workflow

1. Create a branch from `main` (the default branch is protected; every change goes in
   through a pull request).
2. Make your change, with tests.
3. Run what CI runs:
   ```bash
   ruff format --check geomech/ tests/ && ruff check geomech/ tests/
   mypy geomech --ignore-missing-imports
   python -m pytest tests/ -q
   python scripts/regenerate_generated.py --check
   ```
4. Open a pull request against `main`. CI runs lint and the tests on Python 3.10–3.12.

## Rules that keep the equations correct

- **Test the actual equation.** A test of a derived system asserts the full equation
  against a hand derivation or the textbook (Lee, Leok and McClamroch, *Global
  Formulations of Lagrangian and Hamiltonian Dynamics on Manifolds*; cite the equation
  number in the docstring), not a superficial property such as "contains J".
- **Check it numerically.** Every pinned system also passes `geomech.numeric.check_eom`
  in `tests/test_eom_oracle.py`. Add new systems there.
- **Do not edit a passing test to make a change pass.** A failing test means the change
  is wrong, unless the test itself is shown to be wrong; say so in the pull request.
- **Rewrite rules are fuzzed.** After changing anything under `geomech/core/`, run
  `python scripts/engine_fuzz.py --quick`; it must report 0 wrong-value findings.
- **Printed equations.** A change that reprints pinned equation strings goes through
  `scripts/pinned_snapshot.py` (snapshot before, compare after, values must be equal).

`docs/pipeline.md`, `docs/conventions.md` and `docs/troubleshooting.md` explain the
pipeline, the sign conventions and every error message.

## Code style

- `ruff format` and `ruff check` (line length 99). `E741`, `E731` and `E701` are ignored on
  purpose for math variable names and one-line `match` cases.
- Expression nodes are immutable; build a new node instead of changing one.
- Every rewrite goes through `rule(name, before, after)` and is listed with its hypothesis
  in the module docstring (`tests/test_rule_names.py` checks this).

## Reporting bugs

Open an issue with the bug template. The most useful report is a minimal system
definition, what `compute_eom` (or the failing call) returned, and the output of
`check_eom` for it.

## Code of conduct

This project follows the [Contributor Covenant](CODE_OF_CONDUCT.md).
