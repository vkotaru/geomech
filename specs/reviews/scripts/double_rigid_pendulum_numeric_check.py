"""A2: numeric verification of the SO(3)xSO(3) double rigid pendulum against a hand derivation,
plus collect() edge cases for the new (MVMul, MVMul) peel branch."""

import os
import sys

# Check the geomech of the repository this script lives in, not the installed one.
sys.path.insert(0, os.path.abspath(os.path.join(__file__, *[os.pardir] * 4)))

import numpy as np  # noqa: E402

from geomech import (  # noqa: E402
    S2,
    SO3,
    Add,
    Cross,
    Dot,
    Hat,
    MAdd,
    Matrix,
    MMMul,
    Mul,
    MVMul,
    Scalar,
    SMMul,
    SVMul,
    SystemVariables,
    TimeDerivative,
    Transpose,
    VAdd,
    Variation,
    Vector,
    Vee,
    VVMul,
    compute_eom,
    getScalars,
    to_standard_form,
)
from geomech.core.math.collecting import collect

rng = np.random.default_rng(7)


def hat(v):
    v = np.asarray(v).reshape(3)
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], dtype=float)


def vee(M):
    return np.array([[M[2, 1]], [M[0, 2]], [M[1, 0]]])


def rand_rot():
    A = rng.normal(size=(3, 3))
    Q, R = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


class Evaluator:
    def __init__(self):
        self.env = {}

    def leaf(self, node):
        key = str(node)
        if key in self.env:
            return self.env[key]
        if isinstance(node, Scalar):
            if node.value is not None:
                val = float(node.value)
            else:
                val = float(rng.normal())
        elif isinstance(node, SO3):
            val = rand_rot()
        elif isinstance(node, S2):
            v = rng.normal(size=(3, 1))
            val = v / np.linalg.norm(v)
        elif isinstance(node, Matrix):
            if key == "I":
                val = np.eye(3)
            elif key == "0":
                val = np.zeros((3, 3))
            else:
                A = rng.normal(size=(3, 3))
                val = A + A.T if node.is_symmetric else A
        elif isinstance(node, Vector):
            if key in ("0v", "0"):
                val = np.zeros((3, 1))
            else:
                val = rng.normal(size=(3, 1))
        else:
            raise TypeError(f"unknown leaf {type(node)} {key}")
        self.env[key] = val
        return val

    def __call__(self, node):
        if isinstance(node, (Add, VAdd, MAdd)):
            return sum(self(n) for n in node.nodes)
        if isinstance(node, Mul):
            return self(node.left) * self(node.right)
        if isinstance(node, SVMul):
            return self(node.left) * self(node.right)
        if isinstance(node, SMMul):
            return self(node.left) * self(node.right)
        if isinstance(node, (MVMul, MMMul, VVMul)):
            out = self(node.left) @ self(node.right)
            if out.shape == (1, 1):
                return float(out[0, 0])
            return out
        if isinstance(node, Dot):
            return float((self(node.left).T @ self(node.right))[0, 0])
        if isinstance(node, Cross):
            a, b = self(node.left).reshape(3), self(node.right).reshape(3)
            return np.cross(a, b).reshape(3, 1)
        if isinstance(node, Hat):
            return hat(self(node.expr))
        if isinstance(node, Vee):
            return vee(self(node.expr))
        if isinstance(node, Transpose):
            return self(node.expr).T
        if isinstance(node, TimeDerivative):
            # derivative of a leaf: independent random value keyed by its printed name
            inner = node.expr
            if isinstance(inner, (Scalar, Vector, Matrix)):
                key = str(node)
                if key not in self.env:
                    if isinstance(inner, Scalar):
                        self.env[key] = float(rng.normal())
                    elif isinstance(inner, Vector):
                        self.env[key] = rng.normal(size=(3, 1))
                    else:
                        self.env[key] = rng.normal(size=(3, 3))
                return self.env[key]
            raise ValueError(f"TimeDerivative of non-leaf survived: {node}")
        if isinstance(node, Variation):
            raise ValueError(f"Variation survived: {node}")
        return self.leaf(node)


# ---------------------------------------------------------------------------
# Double rigid pendulum (same construction as tests/test_manifold_pipelines.py)
# ---------------------------------------------------------------------------
HALF = Scalar("0.5", value=0.5, attr=["Constant"])
J1 = Matrix("J1", attr=["Constant", "SymmetricMatrix"])
J2 = Matrix("J2", attr=["Constant", "SymmetricMatrix"])
rho1 = Vector("\\rho_1", attr=["Constant"])
rho2 = Vector("\\rho_2", attr=["Constant"])
l1 = Vector("l_1", attr=["Constant"])
m1, m2, g = getScalars("m1 m2 g", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
M1, M2 = Vector("M1"), Vector("M2")
R1, R2 = SO3("R1"), SO3("R2")
Om1, Om2 = R1.get_tangent_vector(), R2.get_tangent_vector()
eta1, eta2 = R1.get_variation_vector(), R2.get_variation_vector()
x1 = R1 * rho1
x2 = R1 * l1 + R2 * rho2
v1, v2 = x1.t_diff(), x2.t_diff()
KE = (
    Dot(Om1, J1 * Om1) * HALF
    + Dot(Om2, J2 * Om2) * HALF
    + Dot(v1, v1) * m1 * HALF
    + Dot(v2, v2) * m2 * HALF
)
PE = m1 * g * Dot(x1, e3) + m2 * g * Dot(x2, e3)
variables = SystemVariables(matrices=[R1, R2])
eom = compute_eom(KE - PE, Dot(eta1, M1) + Dot(eta2, M2), variables)
sf = to_standard_form(eom, variables, [M1, M2])

print("eta1 eqn:", eom[str(eta1)][1])
print("eta2 eqn:", eom[str(eta2)][1])
print("M1:", {k: str(v) for k, v in sf[str(eta1)].M.items()})
print("M2:", {k: str(v) for k, v in sf[str(eta2)].M.items()})
print("f1:", sf[str(eta1)].f)
print("f2:", sf[str(eta2)].f)

worst = 0.0
for trial in range(5):
    ev = Evaluator()
    E1_lib = ev(eom[str(eta1)][1]).reshape(3)
    E2_lib = ev(eom[str(eta2)][1]).reshape(3)
    env = ev.env
    # pull the same random values
    nR1, nR2 = env["R1"], env["R2"]
    nJ1, nJ2 = env["J1"], env["J2"]
    nr1, nr2, nl1 = env["\\rho_1"], env["\\rho_2"], env["l_1"]
    nm1, nm2, ng = env["m1"], env["m2"], env["g"]
    ne3 = env["e3"]
    nM1, nM2 = env["M1"], env["M2"]
    W1, W2 = env[str(Om1)], env[str(Om2)]
    dW1 = env[str(TimeDerivative(Om1))]
    dW2 = env[str(TimeDerivative(Om2))]

    # ---- hand derivation (library convention: coefficient of eta after IBP) ----
    # p_i = dL/dOmega_i
    Jeff1 = nJ1 - nm1 * hat(nr1) @ hat(nr1) - nm2 * hat(nl1) @ hat(nl1)
    Jeff2 = nJ2 - nm2 * hat(nr2) @ hat(nr2)
    C12 = hat(nl1) @ nR1.T @ nR2 @ hat(nr2)  # l1^ R1' R2 rho2^
    C21 = hat(nr2) @ nR2.T @ nR1 @ hat(nl1)  # = C12^T
    p1 = Jeff1 @ W1 - nm2 * C12 @ W2
    p2 = Jeff2 @ W2 - nm2 * C21 @ W1
    # d/dt R1' = -W1^ R1',  d/dt R2 = R2 W2^
    dC12 = hat(nl1) @ (-hat(W1) @ nR1.T @ nR2 @ hat(nr2) + nR1.T @ nR2 @ hat(W2) @ hat(nr2))
    dC21 = hat(nr2) @ (-hat(W2) @ nR2.T @ nR1 @ hat(nl1) + nR2.T @ nR1 @ hat(W1) @ hat(nl1))
    dp1 = Jeff1 @ dW1 - nm2 * (dC12 @ W2 + C12 @ dW2)
    dp2 = Jeff2 @ dW2 - nm2 * (dC21 @ W1 + C21 @ dW1)
    a = np.cross(W1.reshape(3), nl1.reshape(3)).reshape(3, 1)
    b = np.cross(W2.reshape(3), nr2.reshape(3)).reshape(3, 1)

    def cr(u, v):
        return np.cross(u.reshape(3), v.reshape(3)).reshape(3, 1)

    E1 = (
        -dp1
        + cr(p1, W1)
        + nm2 * cr(a, nR1.T @ nR2 @ b)
        - nm1 * ng * cr(nr1, nR1.T @ ne3)
        - nm2 * ng * cr(nl1, nR1.T @ ne3)
        + nM1
    ).reshape(3)
    E2 = (
        -dp2 + cr(p2, W2) + nm2 * cr(b, nR2.T @ nR1 @ a) - nm2 * ng * cr(nr2, nR2.T @ ne3) + nM2
    ).reshape(3)
    err1 = np.max(np.abs(E1 - E1_lib)) / (1 + np.max(np.abs(E1)))
    err2 = np.max(np.abs(E2 - E2_lib)) / (1 + np.max(np.abs(E2)))
    worst = max(worst, err1, err2)
    # M blocks
    dOm1, dOm2 = str(TimeDerivative(Om1)), str(TimeDerivative(Om2))
    eM1 = {k: ev(v) for k, v in sf[str(eta1)].M.items()}
    eM2 = {k: ev(v) for k, v in sf[str(eta2)].M.items()}
    errM = max(
        np.max(np.abs(eM1[dOm1] + Jeff1)),
        np.max(np.abs(eM1[dOm2] - nm2 * C12)),
        np.max(np.abs(eM2[dOm1] - nm2 * C21)),
        np.max(np.abs(eM2[dOm2] + Jeff2)),
        np.max(np.abs(eM1[dOm2] - eM2[dOm1].T)),
    )
    print(f"trial {trial}: rel err eta1={err1:.2e} eta2={err2:.2e}  M-block err={errM:.2e}")
    worst = max(worst, errM)

print("WORST:", worst)
assert worst < 1e-9, "double pendulum does NOT match hand derivation"
print("double rigid pendulum: library == hand derivation (full eta1/eta2 equations and M blocks)")

# ---------------------------------------------------------------------------
# collect() edge cases for the new peel branch
# ---------------------------------------------------------------------------
print("\n--- collect edge cases ---")
sys.setrecursionlimit(400)
vec, av, bv = Vector("v"), Vector("a"), Vector("b")
A, B = Matrix("A"), Matrix("B")
cases = {
    "Dot(A (v x a), B b)": Dot(MVMul(A, Cross(vec, av)), MVMul(B, bv)),
    "Dot(A (a x v), B b)": Dot(MVMul(A, Cross(av, vec)), MVMul(B, bv)),
    "Dot(A (B (v x a)), A b)  (nested, same R outer)": Dot(
        MVMul(A, MVMul(B, Cross(vec, av))), MVMul(A, bv)
    ),
    "Dot(A (B (v x a)), B b)": Dot(MVMul(A, MVMul(B, Cross(vec, av))), MVMul(B, bv)),
    "Dot(A (v s), B b)  (SVMul inside)": Dot(MVMul(A, SVMul(vec, Scalar("s"))), MVMul(B, bv)),
    "Dot(A (v + a), B b)  (VAdd inside)": Dot(MVMul(A, VAdd(vec, av)), MVMul(B, bv)),
    "PRE-EXISTING: Dot(v s, B b)": Dot(SVMul(vec, Scalar("s")), MVMul(B, bv)),
    "both sides have v: Dot(A v, B v)": Dot(MVMul(A, vec), MVMul(B, vec)),
}
for label, e in cases.items():
    try:
        out = collect(e, vec)
        ok = isinstance(out, Dot) and out.left == vec
        print(f"{label:50s} -> {out}   [v on left: {ok}]")
    except RecursionError:
        print(f"{label:50s} -> RecursionError")
    except Exception as ex:  # noqa: BLE001
        print(f"{label:50s} -> {type(ex).__name__}: {ex}")
