"""Models with a symbolic number of bodies, written with ``Sum``.

Each function returns an ``IndexedModel``: the Lagrangian and virtual work as sums over
``i = 1, …, n``, the configuration variables and the inputs. ``compute_indexed_eom``
derives the equation of the i-th body for every n at once; ``model.unroll(size)`` gives
the ordinary Lagrangian, virtual work, variables and inputs for one n, against which
the tests check the generic equations with ``check_eom``. The models are the generic-n
forms of ``n_link_spherical_chain``, ``n_link_rigid_chain`` (with body moments as
inputs), ``quadrotor_flexible_cable`` and the two ``multi_quadrotor`` systems.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from geomech.core.base import S2, SO3, Matrix, Scalar, Vector
from geomech.core.base.constants import HALF
from geomech.core.operations import Cross, Dot
from geomech.indexed.core import Family, IndexedExpr, Sum, indices, n
from geomech.indexed.eom import IndexedVariables, compute_indexed_eom

__all__ = [
    "IndexedModel",
    "spherical_chain",
    "rigid_chain",
    "quadrotor_flexible_cable",
    "multi_quadrotor_point_mass",
    "multi_quadrotor_rigid_load",
    "MODELS",
]


@dataclass
class IndexedModel:
    """A model with families of n bodies: ``L``, ``δW``, variables and inputs."""

    name: str
    configuration: str
    lagrangian: IndexedExpr
    inf_work: IndexedExpr
    variables: IndexedVariables
    input_families: list = field(default_factory=list)
    inputs: list = field(default_factory=list)

    def eom(self):
        """``compute_indexed_eom`` of the model."""
        return compute_indexed_eom(self.lagrangian, self.inf_work, self.variables)

    def unroll(self, size: int):
        """``(L, δW, SystemVariables, inputs)`` of the model with ``n = size``."""
        inputs = list(self.inputs)
        for fam in self.input_families:
            inputs += [fam[c] for c in range(1, size + 1)]
        return (
            _instantiate(self.lagrangian, size),
            _instantiate(self.inf_work, size),
            self.variables.members(size),
            inputs,
        )


def _instantiate(expr, size):
    return expr.instantiate(size) if isinstance(expr, IndexedExpr) else expr


def _const(names):
    return [Scalar(s, attr=["Constant"]) for s in names.split()]


def _inertia_family(name):
    return Family(Matrix, name, attr=["Constant", "SymmetricMatrix"])


def spherical_chain() -> IndexedModel:
    """Chain of n spherical pendulums on (S²)ⁿ (Lee et al. §5.3, eq. 5.14).

    ``x_i = Σ_{j≤i} l_j q_j``; ``L = Σ_i ½ m_i ẋ_i·ẋ_i − Σ_i m_i g e₃·x_i``;
    ``δW = Σ_i ξ_i·M_i``.
    """
    i, j = indices("i j")
    (g,) = _const("g")
    e3 = Vector("e3", attr=["Constant"])
    q = Family(S2, "q")
    m, l = Family(Scalar, "m", ["Constant"]), Family(Scalar, "l", ["Constant"])
    M = Family(Vector, "M")

    def x(k):
        return Sum(j, 1, k, q[j] * l[j])

    def v(k):
        return x(k).t_diff()

    L = Sum(i, 1, n, HALF * m[i] * v(i).dot(v(i))) - Sum(i, 1, n, m[i] * g * x(i).dot(e3))
    dW = Sum(i, 1, n, Dot(q[i].get_variation_vector(), M[i]))
    return IndexedModel("spherical_chain", "(S²)ⁿ", L, dW, IndexedVariables(families=[q]), [M])


def rigid_chain() -> IndexedModel:
    """Chain of n rigid links with spherical joints on SO(3)ⁿ.

    Joint ``p_i = Σ_{j<i} R_j l_j``, centre of mass ``x_i = p_i + R_i ρ_i``;
    ``L = Σ_i ½ Ω_i·J_iΩ_i + Σ_i ½ m_i ẋ_i·ẋ_i − Σ_i m_i g e₃·x_i``; body moments
    ``δW = Σ_i η_i·M_i``.
    """
    i, j = indices("i j")
    (g,) = _const("g")
    e3 = Vector("e3", attr=["Constant"])
    R = Family(SO3, "R")
    J = _inertia_family("J")
    m = Family(Scalar, "m", ["Constant"])
    rho, l = Family(Vector, "rho", ["Constant"]), Family(Vector, "l", ["Constant"])
    M = Family(Vector, "M")

    def x(k):
        return Sum(j, 1, k - 1, R[j] * l[j]) + R[k] * rho[k]

    def omega(k):
        return R[k].get_tangent_vector()

    KE = Sum(i, 1, n, HALF * Dot(omega(i), J[i] * omega(i)))
    KE = KE + Sum(i, 1, n, HALF * m[i] * x(i).t_diff().dot(x(i).t_diff()))
    PE = Sum(i, 1, n, m[i] * g * x(i).dot(e3))
    dW = Sum(i, 1, n, Dot(R[i].get_variation_vector(), M[i]))
    return IndexedModel(
        "rigid_chain", "(SO(3))ⁿ", KE - PE, dW, IndexedVariables(families=[R]), [M]
    )


def quadrotor_flexible_cable() -> IndexedModel:
    """Quadrotor with a load on a cable of n links, on R³ × SO(3) × (S²)ⁿ.

    ``x_i = x_Q + Σ_{j≤i} l_j q_j``; ``L = ½ m_Q ẋ_Q·ẋ_Q + ½ Ω·JΩ − m_Q g e₃·x_Q
    + Σ_i (½ m_i ẋ_i·ẋ_i − m_i g e₃·x_i)``; ``δW = f (R e₃)·δx_Q + η·M``.
    """
    i, j = indices("i j")
    mQ, g, f = Scalar("m_Q", attr=["Constant"]), Scalar("g", attr=["Constant"]), Scalar("f")
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    e3 = Vector("e3", attr=["Constant"])
    xQ, R, M = Vector("x_Q"), SO3("R"), Vector("M")
    q = Family(S2, "q")
    m, l = Family(Scalar, "m", ["Constant"]), Family(Scalar, "l", ["Constant"])
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    vQ = xQ.t_diff()

    def x(k):
        return xQ + Sum(j, 1, k, q[j] * l[j])

    L = HALF * mQ * Dot(vQ, vQ) + HALF * Dot(Om, J * Om) - mQ * g * Dot(e3, xQ)
    L = L + Sum(i, 1, n, HALF * m[i] * x(i).t_diff().dot(x(i).t_diff()))
    L = L - Sum(i, 1, n, m[i] * g * x(i).dot(e3))
    dW = Dot(xQ.delta(), (R * e3) * f) + Dot(eta, M)
    return IndexedModel(
        "quadrotor_flexible_cable",
        "R³ × SO(3) × (S²)ⁿ",
        L,
        dW,
        IndexedVariables(families=[q], vectors=[xQ], matrices=[R]),
        [],
        [f, M],
    )


def _team():
    m, l = Family(Scalar, "m", ["Constant"]), Family(Scalar, "l", ["Constant"])
    J = _inertia_family("J")
    R, q = Family(SO3, "R"), Family(S2, "q")
    f, M = Family(Scalar, "f"), Family(Vector, "M")
    return m, l, J, R, q, f, M


def multi_quadrotor_point_mass() -> IndexedModel:
    """n quadrotors carrying a point-mass load, on R³ × (SO(3) × S²)ⁿ.

    ``x_Qi = x_L − l_i q_i``; ``L = ½ m_L ẋ_L·ẋ_L − m_L g e₃·x_L
    + Σ_i (½ m_i ẋ_Qi·ẋ_Qi + ½ Ω_i·J_iΩ_i − m_i g e₃·x_Qi)``;
    ``δW = Σ_i (f_i (R_i e₃)·(δx_L − l_i ξ_i × q_i) + η_i·M_i)``.
    """
    (i,) = indices("i")
    mL, g = _const("m_L g")
    e3 = Vector("e3", attr=["Constant"])
    xL = Vector("x_L")
    m, l, J, R, q, f, M = _team()
    vL = xL.t_diff()
    xQ = xL - q[i] * l[i]
    Om = R[i].get_tangent_vector()
    L = HALF * mL * Dot(vL, vL) - mL * g * Dot(e3, xL)
    L = L + Sum(
        i,
        1,
        n,
        HALF * m[i] * Dot(xQ.t_diff(), xQ.t_diff())
        + HALF * Dot(Om, J[i] * Om)
        - m[i] * g * Dot(e3, xQ),
    )
    dxQ = xL.delta() - Cross(q[i].get_variation_vector(), q[i]) * l[i]
    dW = Sum(i, 1, n, Dot(dxQ, (R[i] * e3) * f[i]) + Dot(R[i].get_variation_vector(), M[i]))
    return IndexedModel(
        "multi_quadrotor_point_mass",
        "R³ × (SO(3) × S²)ⁿ",
        L,
        dW,
        IndexedVariables(families=[q, R], vectors=[xL]),
        [f, M],
    )


def multi_quadrotor_rigid_load() -> IndexedModel:
    """n quadrotors carrying a rigid-body load, on R³ × SO(3) × (SO(3) × S²)ⁿ.

    Cable i is attached at ``x_L + R_L ρ_i``; ``x_Qi = x_L + R_L ρ_i − l_i q_i``.
    """
    (i,) = indices("i")
    mL, g = _const("m_L g")
    JL = Matrix("J_L", attr=["Constant", "SymmetricMatrix"])
    e3 = Vector("e3", attr=["Constant"])
    xL, RL = Vector("x_L"), SO3("R_L")
    OmL, etaL = RL.get_tangent_vector(), RL.get_variation_vector()
    rho = Family(Vector, "rho", ["Constant"])
    m, l, J, R, q, f, M = _team()
    vL = xL.t_diff()
    xQ = xL + RL * rho[i] - q[i] * l[i]
    Om = R[i].get_tangent_vector()
    L = HALF * mL * Dot(vL, vL) + HALF * Dot(OmL, JL * OmL) - mL * g * Dot(e3, xL)
    L = L + Sum(
        i,
        1,
        n,
        HALF * m[i] * Dot(xQ.t_diff(), xQ.t_diff())
        + HALF * Dot(Om, J[i] * Om)
        - m[i] * g * Dot(e3, xQ),
    )
    dxQ = xL.delta() + RL * Cross(etaL, rho[i]) - Cross(q[i].get_variation_vector(), q[i]) * l[i]
    dW = Sum(i, 1, n, Dot(dxQ, (R[i] * e3) * f[i]) + Dot(R[i].get_variation_vector(), M[i]))
    return IndexedModel(
        "multi_quadrotor_rigid_load",
        "R³ × SO(3) × (SO(3) × S²)ⁿ",
        L,
        dW,
        IndexedVariables(families=[q, R], vectors=[xL], matrices=[RL]),
        [f, M],
    )


MODELS = {
    fn.__name__: fn
    for fn in (
        spherical_chain,
        rigid_chain,
        quadrotor_flexible_cable,
        multi_quadrotor_point_mass,
        multi_quadrotor_rigid_load,
    )
}
