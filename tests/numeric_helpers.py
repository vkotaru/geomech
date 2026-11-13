"""Re-export of ``geomech.numeric`` kept so that the test suite reads as before.

The evaluator and the finite-difference helpers now live in the package
(``geomech/numeric.py``) so that they can be imported from anywhere; this module only
forwards the names the tests use.
"""

from __future__ import annotations

from geomech.numeric import (
    Evaluator,
    assert_eom,
    check_eom,
    evaluate,
    expm_so3,
    finite_difference_check,
    hat,
    linearized_rhs,
    perturbed_state,
    rand_rot,
    sample_direction,
    sample_state,
    vee,
)

__all__ = [
    "Evaluator",
    "assert_eom",
    "check_eom",
    "evaluate",
    "expm_so3",
    "finite_difference_check",
    "hat",
    "linearized_rhs",
    "perturbed_state",
    "rand_rot",
    "sample_direction",
    "sample_state",
    "vee",
]
