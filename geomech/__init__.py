"""
geomech is a python library for symbolic computation of Dynamics on Manifolds.
"""

import sys

if sys.version_info < (3, 10):
    raise ImportError("Python version >= 3.10 is required")
del sys

from importlib.metadata import PackageNotFoundError as _PackageNotFoundError
from importlib.metadata import version as _version

try:
    __version__ = _version("geomech")
except _PackageNotFoundError:  # running from a source tree that is not installed
    __version__ = "unknown"
del _version, _PackageNotFoundError

# --- Layer 1: types & leaf nodes ---
from geomech.core.base import (
    S2,
    SO3,
    TS2,
    TSO3,
    Expr,
    ExprFlags,
    ExprType,
    I,
    IdentityMatrix,
    ManifoldInfo,
    Matrix,
    MatrixExpr,
    Number,
    O,
    One,
    Scalar,
    ScalarExpr,
    SkewSymmMatrix,
    Vector,
    VectorExpr,
    Zero,
    ZeroMatrix,
    ZeroVector,
    getMatrices,
    getScalars,
    getVectors,
)

# --- Layer 4: symbolic math ---
from geomech.core.math.collecting import collect
from geomech.core.math.extract import extract_linear_coeff
from geomech.core.math.ibp import integrate_by_parts

# --- Layer 2: expression tree nodes ---
from geomech.core.operations import (
    Add,
    Cross,
    Dot,
    Hat,
    MAdd,
    MMMul,
    Mul,
    MVMul,
    SMMul,
    SVMul,
    TimeDerivative,
    TimeIntegral,
    Transpose,
    VAdd,
    Variation,
    Vee,
    VVMul,
)

# --- Layer 3: algorithms ---
from geomech.core.transformations.expansion import expand
from geomech.core.transformations.manifold_rules import apply_manifold_rules
from geomech.core.transformations.pulling import pull
from geomech.core.transformations.simplification import full_simplify, simplify
from geomech.core.transformations.vector_identities import vector_rules

# --- Layer 5: dynamics pipeline ---
from geomech.dynamics import (
    LinearizedSystem,
    StandardFormEquation,
    SystemVariables,
    compute_eom,
    linearize,
    separate_variations,
    to_standard_form,
)

# --- Numeric oracle (imports only layers 1-2 itself, so it may sit anywhere after them) ---
from geomech.numeric import (
    Evaluator,
    assert_eom,
    check_eom,
    evaluate,
    finite_difference_check,
)

# --- Utilities ---
from geomech.utils import (
    display_eom,
    display_latex,
    display_standard_form,
    eom_to_latex,
    print_eom,
    print_tree,
    render_eom,
    to_latex,
    tree_str,
)
from geomech.utils.errors import (
    AlgebraicError,
    ExpressionMismatchError,
    InvariantError,
    NonLinearError,
    SizeMismatchError,
)
from geomech.utils.invariants import check_tree
from geomech.utils.trace import EOMTrace, LinearizationTrace, set_trace
from geomech.utils.trace import configure_from_env as _configure_from_env

# GEOMECH_DEBUG=1: the "geomech" logger at DEBUG and check_tree at pipeline entry/exit;
# GEOMECH_TRACE=1: log every named rewrite rule (see geomech/utils/trace.py)
_configure_from_env()
