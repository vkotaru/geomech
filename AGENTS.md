# geomech: notes for coding agents

`geomech` derives coordinate-free equations of motion on R, R³, S¹, S², SO(3) and their
products from a Lagrangian and virtual work. Documentation: https://geomech.readthedocs.io

## Deriving equations of motion

To derive, check, linearize, generate code for or simulate a mechanical system, follow the
skill in [`skills/geomech/SKILL.md`](skills/geomech/SKILL.md) (API and error messages in
[`skills/geomech/reference.md`](skills/geomech/reference.md)). Start from
`geomech.systems.get(name)` when the system is in the library, and always run `check_eom`
and report its error. The skill is the only copy of these instructions; to use it in
another project, copy `skills/geomech/` into the directory your agent reads (see the
README, "Using geomech with AI coding agents").

## Working on the repository

Setup, the checks CI runs, and the rules for new systems are in
[`CONTRIBUTING.md`](CONTRIBUTING.md). In short:

```bash
python -m pip install -e ".[dev]"
ruff format --check geomech/ tests/ && ruff check geomech/ tests/
mypy geomech --ignore-missing-imports
python -m pytest tests/ -q
```

- Tests check the derived equations against the literature or a numerical model. A failing
  test means the change is wrong; do not edit a passing test to make a change pass.
