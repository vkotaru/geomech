"""Unit tests of the validation layer inside ``linearize``.

End-to-end behaviour is covered in ``tests/test_vbl.py`` and
``tests/test_vbl_adversarial.py``; this file covers the internal pieces those tests cannot
reach on their own:

- ``_partition_terms`` rejects a term carrying two perturbation symbols even when
  ``expand`` could not split it (with the ``VVMul`` expand rule in place the only way to
  reach the check is to call it directly);
- ``_occurs_outside`` skips the owner subtree so ``δω`` inside the owner ``d/dt(δω)`` is
  not a second symbol;
- blocks are routed to M / A / B by position, and the resulting ``LinearizedSystem`` is
  dense over ``inputs`` and ``states``.
"""

import pytest

from geomech import (
    S2,
    SO3,
    Cross,
    Dot,
    Matrix,
    MVMul,
    SystemVariables,
    TimeDerivative,
    Transpose,
    Variation,
    Vector,
    VVMul,
    compute_eom,
    linearize,
)
from geomech.dynamics.vbl import _occurs_outside, _partition_terms
from geomech.utils.errors import AlgebraicError
from tests import vbl_systems

q = S2("q")
om, xi = q.get_tangent_vector(), q.get_variation_vector()
dom = Variation(om)
ddom = TimeDerivative(dom)
e3 = Vector("e3", attr=["Constant"])


class TestPartitionExactlyOne:
    def test_term_with_two_symbols_raises(self):
        """A single undistributed term (ξ × q + δω) e₃ᵀ e₃ must not be bucketed."""
        term = MVMul(VVMul(Cross(xi, q) + dom, Transpose(e3)), e3)
        with pytest.raises(AlgebraicError, match="carries both"):
            _partition_terms(term, [ddom, dom, xi], "row")

    def test_nonlinear_term_raises(self):
        """δω (ξ·e₃) carries two symbols even though every factor is distributed."""
        term = dom * Dot(xi, e3)
        with pytest.raises(AlgebraicError, match="carries both"):
            _partition_terms(term, [ddom, dom, xi], "row")

    def test_owner_inside_calculus_node_is_not_a_second_symbol(self):
        """J d/dt(δω) is owned by d/dt(δω); the δω inside it must not trip the check."""
        J = Matrix("J", attr=["Constant", "Symmetric"])
        term = MVMul(J, ddom)
        buckets = _partition_terms(term, [ddom, dom, xi], "row")
        assert [len(b) for b in buckets] == [1, 0, 0]

    def test_terms_land_in_their_own_buckets(self):
        dE = MVMul(Matrix("J"), ddom) + Cross(dom, e3) + Cross(xi, q)
        buckets = _partition_terms(dE, [ddom, dom, xi], "row")
        assert [len(b) for b in buckets] == [1, 1, 1]


class TestOccursOutside:
    def test_owner_subtree_is_skipped(self):
        assert not _occurs_outside(MVMul(Matrix("J"), ddom), dom, ddom)

    def test_symbol_outside_owner_is_found(self):
        assert _occurs_outside(Cross(xi, q) + dom, xi, dom)

    def test_symbol_inside_foreign_calculus_node_counts(self):
        """ξ̇ is not the owner δω and hides ξ: it counts as a second symbol."""
        assert _occurs_outside(Cross(TimeDerivative(xi), q) + dom, xi, dom)


class TestRoutingByPosition:
    def test_blocks_are_dense_over_states_and_inputs(self):
        s = vbl_systems.free_rigid_body()
        eom = compute_eom(s.L, s.dW, s.variables)
        lin = linearize(eom, s.variables, [s.symbols["M"]])
        for row in lin.rows:
            assert set(lin.M[row]) == {str(d) for d in lin.dstates}
            assert set(lin.A[row]) == {str(x) for x in lin.states}
            assert set(lin.B[row]) == {str(u) for u in lin.inputs}
            assert len(lin.B[row]) == len(lin.inputs)

    def test_textbook_block_stays_in_A(self):
        """Hat(JΩ) − Hat(Ω) J is an A block, never a B block, whatever the inputs look like."""
        s = vbl_systems.free_rigid_body()
        eom = compute_eom(s.L, s.dW, s.variables)
        R = s.symbols["R"]
        eta, Om = R.get_variation_vector(), R.get_tangent_vector()
        lin = linearize(eom, s.variables, [s.symbols["M"]])
        assert not lin.A[str(eta)][str(Variation(Om))].is_zero
        assert set(lin.B[str(eta)]) == {str(Variation(s.symbols["M"]))}


class TestInputValidationMessages:
    def test_tangent_vector_input_names_the_problem(self):
        R = SO3("R")
        s = vbl_systems.free_rigid_body()
        eom = compute_eom(s.L, s.dW, s.variables)
        with pytest.raises(ValueError, match="tangent or variation vector"):
            linearize(eom, s.variables, [s.symbols["R"].get_tangent_vector()])
        del R

    def test_input_whose_variation_is_a_state_is_rejected(self):
        """A plain Vector printing like the tangent symbol Ω collides with the state δΩ."""
        s = vbl_systems.free_rigid_body()
        eom = compute_eom(s.L, s.dW, s.variables)
        look_alike = Vector(str(s.symbols["R"].get_tangent_vector()))
        with pytest.raises(ValueError, match="already a perturbation state"):
            linearize(eom, s.variables, [s.symbols["M"], look_alike])

    def test_duplicate_input_message(self):
        s = vbl_systems.free_rigid_body()
        eom = compute_eom(s.L, s.dW, s.variables)
        with pytest.raises(ValueError, match="more than once"):
            linearize(eom, s.variables, [s.symbols["M"], s.symbols["M"]])

    def test_matrix_equation_is_rejected(self):
        R = SO3("R")
        eta = R.get_variation_vector()
        with pytest.raises(TypeError, match="scalar or vector"):
            linearize({str(eta): (eta, R)}, SystemVariables(matrices=[R]))

    def test_type_mismatch_message_names_both_types(self):
        R = SO3("R")
        Om, eta = R.get_tangent_vector(), R.get_variation_vector()
        with pytest.raises(TypeError, match="SCALAR.*VECTOR"):
            linearize({str(eta): (eta, Dot(Om, Om))}, SystemVariables(matrices=[R]))
