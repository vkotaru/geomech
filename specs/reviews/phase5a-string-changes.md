# Phase 5a: changed pinned strings

Every pinned string that phase 5a changed in an existing test, with the numeric check that
licensed the change. Base: `fix/phase4-debuggability` at 47cce2d.

## How the strings were checked

`specs/reviews/scripts/pinned_snapshot.py snapshot` was run on a checkout of 47cce2d. It
recorded every `compute_eom`, `to_standard_form`, `linearize` and `integrate_by_parts`
call made by the test suite (442 distinct calls: pickled arguments plus the ids of the
tests that made them) and a `compute_eom → to_standard_form → linearize` chain for each of
the 67 builder systems of `tests/vbl_systems.py`, `tests/phase1_systems.py` and
`tests/test_eom_oracle.py`; 509 records in all. For each output it stored `str` and the
value at five fixed states (seeds 1001–1005; configuration on the manifold, every other
leaf drawn from a generator seeded by its printed name, so the value does not depend on
how the tree is written), each `linearize` row `−M δẋ + A δx + B δu` along an admissible
direction at the same states, and `check_eom` for each `compute_eom` call.

`pinned_snapshot.py compare` replays the same calls with the phase-5a tree:

| | count |
|---|---|
| outputs whose string is unchanged | 3055 |
| outputs whose string changed | 655 |
| largest relative value difference among those, `max abs(new − old) / max(1, max abs(old))` | 7.3e-16 (649 outputs) |
| `linearize` blocks of three S² systems that differ as matrices | 6 (see below) |
| linearization rows that differ along admissible directions | 0 |
| expected errors still raised | 152 |
| expected errors that now succeed | 2 (see below) |
| `check_eom`: compared / worse than before (above max(old, 1e-6)) | 197 / 0 (largest ratio new/old 1.000000006) |

The six S² blocks are `A[ξ][ξ]` and `A[ξ][δω]` of `s2_anisotropic_qdot`,
`adv_omega_and_qdot_mixed` and `adv_omega_dot_S_qdot` (builder systems; no test pins
them). Their equations `E` print shorter (a tangency identity now fires below an `MVMul`)
and evaluate to the same value on the sphere. The variations of two forms that agree on
the sphere differ by terms proportional to `δq·ω + q·δω`, which vanish along every
admissible perturbation (`ξ ⊥ q`, `q·δω = −ω·(ξ × q)`). Measured: every linearization row
along admissible directions agrees to 1e-15, and the `δω` column difference is of the form
`u qᵀ` (`d (I − qqᵀ)` at most 1.6e-15). The blocks were never unique on S² (the old engine
already applied the same identities at the top level of a term).

The two calls that used to raise are `linearize` on `MVMul(Transpose(Hat(q) + Hat(ω)), ω)`
and `Vee(Hat(q + ω))` (`tests/test_vbl_adversarial.py::TestErrorPathsHeld`):
`(Hat(a) + Hat(b))ᵀ` now simplifies to `−(Hat(a) + Hat(b))` and `Vee` distributes over the
sum, so both linearize. The parametrized test was split; the two cases now assert
`finite_difference_check < 1e-6` (they give 1.5e-10 and 6.6e-11).

Strings in the tests were then updated mechanically from the new output, only where the
old string was recorded by the snapshot and the new one is the replayed output of the same
call; every such test passes with the new string, so each new string is exactly a replayed,
numerically checked output.

## The systems from Lee et al.

| system | equation | old → new length | terms | max rel. diff | new form vs book |
|---|---|---|---|---|---|
| free rigid body (Eq. 6.16) | η | 69 → 67 | 3 → 3 | 0 | closer: `(-1.0)` → `(-1)` |
| rigid pendulum (pivot) | η | 94 → 92 | 4 → 4 | 0 | closer: `(-1.0)` → `(-1)` |
| heavy top (pivot, offset ρ; derived ad hoc, not in the suite) | η | 93 → 91 | 4 → 4 | 0 | closer: `(-1.0)` → `(-1)` |
| spherical pendulum | ξ | 54 → 54 | 3 → 3 | 0 | equal: `ω̇(-1)mll` → `ω̇mll(-1)` (the number last, as in every other term) |
| double spherical pendulum (derived ad hoc, not in the suite) | ξ₁, ξ₂ | 683 → 445, 605 → 367 | 14 → 10, 12 → 8 | 4.4e-16 | closer: the `m₂l₁l₂` and `m₂l₂l₁` copies of each coupling term merge, `(-0.5) + (-0.5) → (-1)` |
| quadrotor | x, η | 50 → 50, 69 → 67 | 3, 3 | 0 | equal / closer |
| quadrotor with load | x_L, ξ, η | 180 → 180, 285 → 285, 69 → 67 | 7, 7, 3 | 1.2e-16 | equal / equal / closer |
| double rigid pendulum | η₁, η₂ | 700 → 698, 561 → 559 | 14, 11 | 0 | closer: `(-1.0)` → `(-1)` |
| elastic spherical pendulum | l, ξ | 98 → 96, 88 → 88 | 6, 4 | 0 | closer (`kl(-1.0)` → `kl(-1)`) / equal |
| rotating sphere (Eq. 5.43), ω- and q̇-form | ξ | 319 → 317, 333 → 331 | 8, 8 | 8.0e-17 | closer: `(2.0)` → `(2)`; same terms |

"Equal" means the same terms with the factors in a different order (numbers last); no
system got longer. The double spherical pendulum and the heavy top were derived with both
trees (`check_eom` 4.37e-8 and 3.7e-9 before and after). `check_eom` of every other system
above is unchanged to the printed digits.

## Changed pinned strings

One row per changed hunk of an existing test that touches a string literal (pinned
strings, and the docstring lines added where the printed form now differs from the hand
derivation). "max rel. diff" is the largest relative value difference over every pipeline
output of the test (0 when the outputs only changed form); "same output in a builder
chain" marks a call made in `setup_method`, whose identical string is an output of a
builder-system chain; "n/a" is a `simplify` unit test (`Mul(Mul(2, 3), a)`), not a
pipeline call. Structural asserts that moved with the numbers-last order
(`tests/test_ibp.py`: `result.left.value == -1` → `result.right.value == -1`;
`tests/test_simplify.py::test_chain_mul_numerics`), the D1 assert in
`tests/test_phase1_adversarial.py::TestOracle` (`"m" not in` → `"m" in`) and the three
strict xfails of `tests/test_scalar_adversarial.py` that now pass (markers removed:
`test_numeric_coefficients_have_one_type`, `test_f_combines_like_terms`,
`test_xi_f_exact_structure`) are listed in the implementation report.

| test | old | new | max rel. diff |
|---|---|---|---|
| test_collect_fix.py::TestGyroscopicCouplingSO3 | `"(J\\frac{d}{dt}(\\Omega_{R})(-1.0)+Cross(\\Omega_{R},J\\Omega_{R})(-1)"` | `"(J\\frac{d}{dt}(\\Omega_{R})(-1)+Cross(\\Omega_{R},J\\Omega_{R})(-1)"` | 0.0e+00 |
| test_collect_fix.py::TestInertialAngularMomentumSO3 | `"(J\\frac{d}{dt}(\\Omega_{R})(-1.0)+Cross(\\Omega_{R},J\\Omega_{R})(-1)""+Cross(J\\Omega_{R},(R)'e3)+JCross(\\Omega_{R},(R)'e3)+Cross(\\Omega_{R},J(R)'e3)(-1)+M)"` | `"(J\\frac{d}{dt}(\\Omega_{R})(-1)+Cross(\\Omega_{R},J\\Omega_{R})(-1)""+Cross(J\\Omega_{R},(R)'e3)+JCross(\\Omega_{R},(R)'e3)+Cross(\\Omega_{R},J(R)'e3)(-1)""+M)"` | 0.0e+00 |
| test_collect_fix.py::TestAnisotropicInertiaQdotForm::test_pinnedequation | `"(Cross(Cross(\\omega_{q},q),ACross(\\omega_{q},q))(-1.0)""+Cross(q,ACross(\\frac{d}{dt}(\\omega_{q}),q))(-1.0)""+Cross(q,ACross(\\omega_{q},Cross(\\omega_{q},q)))(-1.0)""+Cross(\\omega_{q},Cross(q,ACross(\\omega_{q},q)))(2.0)""+Cross(q,Cross(\\omega_{q},ACross(\\omega_{q},q)))(-1)""+Cross(q,e3)mg(-1)+f)"` | `"(Cross(Cross(\\omega_{q},q),ACross(\\omega_{q},q))(-1)""+Cross(q,ACross(\\frac{d}{dt}(\\omega_{q}),q))(-1)""+Cross(q,Aq)Dot(\\omega_{q},\\omega_{q})""+Cross(\\omega_{q},Cross(q,ACross(\\omega_{q},q)))(2)""+Cross(q,Cross(\\omega_{q},ACross(\\omega_{q},q)))(-1)+Cross(q,e3)mg(-1)+f)"` | 3.8e-16 |
| test_collect_fix.py::TestRelativeRateCouplingSO3xSO3::test_pinned_equations | `"(J1\\frac{d}{dt}(\\Omega_{R1})(-1.0)+Cross(\\Omega_{R1},J1\\Omega_{R1})(-1)""+Cross(\\Omega_{R1},(R1)'R2\\Omega_{R2})k""+(R1)'R2Cross(\\Omega_{R2},\\Omega_{R2})(-1)k""+(R1)'R2\\frac{d}{dt}(\\Omega_{R2})(-1)k""+Cross(\\Omega_{R1},(R1)'R2\\Omega_{R2})(-1)k""+Cross((R1)'R2\\Omega_{R2},\\Omega_{R1})(-1)k+M1)"` | `"(J1\\frac{d}{dt}(\\Omega_{R1})(-1)+Cross(\\Omega_{R1},J1\\Omega_{R1})(-1)""+(R1)'R2\\frac{d}{dt}(\\Omega_{R2})k(-1)""+Cross((R1)'R2\\Omega_{R2},\\Omega_{R1})k(-1)+M1)"` | 3.0e-16 |
| test_collect_fix.py::TestRelativeRateCouplingSO3xSO3::test_pinned_equations | `"(J2\\frac{d}{dt}(\\Omega_{R2})(-1.0)+Cross(\\Omega_{R2},J2\\Omega_{R2})(-1)""+Cross(\\Omega_{R2},(R2)'R1\\Omega_{R1})k(2)""+(R2)'R1Cross(\\Omega_{R1},\\Omega_{R1})(-1)k""+(R2)'R1\\frac{d}{dt}(\\Omega_{R1})(-1)k""+Cross(\\Omega_{R2},(R2)'R1\\Omega_{R1})(-1)k+M2)"` | `"(J2\\frac{d}{dt}(\\Omega_{R2})(-1)+Cross(\\Omega_{R2},J2\\Omega_{R2})(-1)""+Cross(\\Omega_{R2},(R2)'R1\\Omega_{R1})k+(R2)'R1\\frac{d}{dt}(\\Omega_{R1})k(-1)""+M2)"` | 3.0e-16 |
| test_collect_fix.py::TestOffsetMassMatrixInertia::test_pinnedequation | `"(Cross(\\rho,JCross(\\frac{d}{dt}(\\Omega_{R}),\\rho))(-1.0)"` | `"(Cross(\\rho,JCross(\\frac{d}{dt}(\\Omega_{R}),\\rho))(-1)"` | 0.0e+00 |
| test_error_context.py::test_row_kinetic_energy_derives | `"+Kx(-1.0)+f)"` | `"+Kx(-1)+f)"` | 0.0e+00 |
| test_ibp.py::TestIntegrateByParts::test_simple_dot | `assert result.left.value == -1assert isinstance(result.right, Dot)assert str(result.right.left) == "eta"assert isinstance(result.right.right, TimeDerivative)` | `assert result.right.value == -1assert isinstance(result.left, Dot)assert str(result.left.left) == "eta"assert isinstance(result.left.right, TimeDerivative)` | 0.0e+00 |
| test_ibp_prefactor.py::TestConstantPrefactorUnchanged::test_full_pipeline | `assert str(out) == f"c(-1)Dot(xi,{DOM})"` | `assert str(out) == f"cDot(xi,{DOM})(-1)"` | 0.0e+00 |
| test_ibp_prefactor.py::TestNonConstantPrefactor::test_single_factor | `assert str(out) == f"(a(-1)Dot(xi,{DOM})+{DA}Dot(xi,omega)(-1))"` | `assert str(out) == f"(aDot(xi,{DOM})(-1)+{DA}Dot(xi,omega)(-1))"` | 0.0e+00 |
| test_ibp_prefactor.py::TestNonConstantPrefactor::test_prefactor_inside_dot_is_pulled_first | `assert str(out) == f"(a(-1)Dot(xi,{DOM})+{DA}Dot(xi,omega)(-1))"` | `assert str(out) == f"(aDot(xi,{DOM})(-1)+{DA}Dot(xi,omega)(-1))"` | 0.0e+00 |
| test_ibp_prefactor.py::TestNonConstantPrefactor::test_two_factors_product_rule | `assert str(out) == (f"(ab(-1)Dot(xi,{DOM})+{DA}bDot(xi,omega)(-1)+a{DB}Dot(xi,omega)(-1))")` | `assert str(out) == (f"(abDot(xi,{DOM})(-1)+{DA}bDot(xi,omega)(-1)+a{DB}Dot(xi,omega)(-1))")` | 0.0e+00 |
| test_ibp_prefactor.py::TestNonConstantPrefactor::test_repeated_factor_is_one_term | `assert str(out) == f"(ll(-1)Dot(xi,{DOM})+l{DL}Dot(xi,omega)(-2))"` | `assert str(out) == f"(llDot(xi,{DOM})(-1)+l{DL}Dot(xi,omega)(-2))"` | 0.0e+00 |
| test_ibp_prefactor.py::TestNonConstantPrefactor::test_constant_rhs_keeps_only_prefactor_rate | `assert str(out) == f"(-1){DA}Dot(xi,e)"` | `assert str(out) == f"{DA}Dot(xi,e)(-1)"` | 0.0e+00 |
| test_ibp_prefactor.py::TestNonConstantPrefactor::test_dot_in_middle_of_chain_not_double_counted | `assert str(out) == (f"(a(-1)Dot(xi,{DOM})b+{DA}bDot(xi,omega)(-1)+a{DB}Dot(xi,omega)(-1))")` | `assert str(out) == (f"(aDot(xi,{DOM})b(-1)+{DA}bDot(xi,omega)(-1)+a{DB}Dot(xi,omega)(-1))")` | 0.0e+00 |
| test_ibp_prefactor.py::TestNonConstantPrefactor::test_mixed_sum | `f"(c(-1)Dot(xi,{DOM})+a(-1)Dot(xi,\\frac{{d}}{{dt}}(w))+{DA}Dot(xi,w)(-1))"` | `f"(cDot(xi,{DOM})(-1)+aDot(xi,\\frac{{d}}{{dt}}(w))(-1)+{DA}Dot(xi,w)(-1))"` | 0.0e+00 |
| test_manifold_pipelines.py::TestSphericalPendulumPipeline::test_eom_str | `expected = "\\frac{d}{dt}(\\omega_{q})(-1)mll+Cross(q,e3)mgl(-1)+f"` | `expected = "\\frac{d}{dt}(\\omega_{q})mll(-1)+Cross(q,e3)mgl(-1)+f"` | 0.0e+00 |
| test_manifold_pipelines.py::TestSphericalPendulumPipeline::test_M | `assert str(self.eq.M[ddw]) == "I(-1)mll"` | `assert str(self.eq.M[ddw]) == "Imll(-1)"` | 0.0e+00 |
| test_manifold_pipelines.py::TestSphericalPendulumQdotPipeline::test_eom_str | `expected = "Cross(q,Cross(\\frac{d}{dt}(\\omega_{q}),q))(-1)llm+Cross(q,e3)mgl(-1)+f"` | `expected = "Cross(q,Cross(\\frac{d}{dt}(\\omega_{q}),q))mll(-1)+Cross(q,e3)mgl(-1)+f"` | 1.4e-16 |
| test_manifold_pipelines.py::TestSphericalPendulumQdotPipeline::test_M | `assert str(self.eq.M[ddw]) == "Hat(q)Hat(q)llm"` | `assert str(self.eq.M[ddw]) == "Hat(q)Hat(q)mll"` | 1.4e-16 |
| test_manifold_pipelines.py::TestFreeRigidBodyPipeline::test_eom_str | `expected = "J\\frac{d}{dt}(\\Omega_{R})(-1.0)+Cross(\\Omega_{R},J\\Omega_{R})(-1)+M"` | `expected = "J\\frac{d}{dt}(\\Omega_{R})(-1)+Cross(\\Omega_{R},J\\Omega_{R})(-1)+M"` | 0.0e+00 |
| test_manifold_pipelines.py::TestFreeRigidBodyPipeline::test_M | `assert str(self.eq.M[ddOm]) == "J(-1.0)"` | `assert str(self.eq.M[ddOm]) == "J(-1)"` | 0.0e+00 |
| test_manifold_pipelines.py::TestRigidPendulumPivotPipeline::test_eom_str | `"J\\frac{d}{dt}(\\Omega_{R})(-1.0)"` | `"J\\frac{d}{dt}(\\Omega_{R})(-1)"` | 0.0e+00 |
| test_manifold_pipelines.py::TestRigidPendulumPivotPipeline::test_M | `assert str(self.eq.M[ddOm]) == "J(-1.0)"` | `assert str(self.eq.M[ddOm]) == "J(-1)"` | 0.0e+00 |
| test_manifold_pipelines.py::TestRigidPendulumCOMPipeline::test_eom_str | `"J\\frac{d}{dt}(\\Omega_{R})(-1.0)"` | `"J\\frac{d}{dt}(\\Omega_{R})(-1)"` | 0.0e+00 |
| test_manifold_pipelines.py::TestRigidPendulumCOMPipeline::test_eom_str | `"+Cross(\\rho,Cross(\\frac{d}{dt}(\\Omega_{R}),\\rho))(-1)m""+Cross(\\Omega_{R},Cross(\\rho,Cross(\\Omega_{R},\\rho)))(-1)m"` | `"+Cross(\\rho,Cross(\\frac{d}{dt}(\\Omega_{R}),\\rho))m(-1)""+Cross(\\Omega_{R},Cross(\\rho,Cross(\\Omega_{R},\\rho)))m(-1)"` | 0.0e+00 |
| test_manifold_pipelines.py::TestRigidPendulumCOMPipeline::test_M | `assert str(self.eq.M[ddOm]) == "(J(-1.0)+Hat(\\rho)Hat(\\rho)m)"` | `assert str(self.eq.M[ddOm]) == "(J(-1)+Hat(\\rho)Hat(\\rho)m)"` | 0.0e+00 |
| test_manifold_pipelines.py::TestMixedR3SO3Pipeline::test_translational_eom | `assert str(eq.M[ddx]) == "I(-1)m"` | `assert str(eq.M[ddx]) == "Im(-1)"` | 0.0e+00 |
| test_manifold_pipelines.py::TestSpringMassPipeline::test_eom_str | `expected = "m\\frac{d}{dt}(\\frac{d}{dt}(x))(-1.0)+kx(-1.0)+f"` | `expected = "m\\frac{d}{dt}(\\frac{d}{dt}(x))(-1)+kx(-1)+f"` | 0.0e+00 |
| test_manifold_pipelines.py::TestSpringMassDamperPipeline::test_eom_str | `expected = "m\\frac{d}{dt}(\\frac{d}{dt}(x))(-1.0)+kx(-1.0)+f+c\\frac{d}{dt}(x)(-1)"` | `expected = "m\\frac{d}{dt}(\\frac{d}{dt}(x))(-1)+kx(-1)+f+c\\frac{d}{dt}(x)(-1)"` | 0.0e+00 |
| test_manifold_pipelines.py::TestDoubleSpringMassPipeline::test_x1_eom_str | `expected = "m1\\frac{d}{dt}(\\frac{d}{dt}(x1))(-1.0)+k1x1(-1.0)+k2x2+k2x1(-1.0)+f1"` | `expected = "m1\\frac{d}{dt}(\\frac{d}{dt}(x1))(-1)+k1x1(-1)+k2x2+k2x1(-1)+f1"` | 0.0e+00 |
| test_manifold_pipelines.py::TestDoubleSpringMassPipeline::test_x2_eom_str | `expected = "m2\\frac{d}{dt}(\\frac{d}{dt}(x2))(-1.0)+k2x2(-1.0)+k2x1+f2"` | `expected = "m2\\frac{d}{dt}(\\frac{d}{dt}(x2))(-1)+k2x2(-1)+k2x1+f2"` | 0.0e+00 |
| test_manifold_pipelines.py::TestDoubleRigidPendulumMassMatrixExact::test_eta1_blocks_exact | `assert str(M[self.dOm1]) == "(J1(-1.0)+Hat(\\rho_1)Hat(\\rho_1)m1+Hat(l_1)Hat(l_1)m2)"` | `assert str(M[self.dOm1]) == "(J1(-1)+Hat(\\rho_1)Hat(\\rho_1)m1+Hat(l_1)Hat(l_1)m2)"` | 0.0e+00 |
| test_manifold_pipelines.py::TestDoubleRigidPendulumMassMatrixExact::test_eta2_blocks_exact | `assert str(M[self.dOm2]) == "(J2(-1.0)+Hat(\\rho_2)Hat(\\rho_2)m2)"` | `assert str(M[self.dOm2]) == "(J2(-1)+Hat(\\rho_2)Hat(\\rho_2)m2)"` | 0.0e+00 |
| test_phase1_adversarial.py::TestOracle::test_scalar_with_value_is_folded_by_engine_and_used_by_evaluator | `` | `"""Decision D1 (phase 5a): a named scalar with a value stays symbolic in theengine (it used to be folded to the literal 2); the evaluator still uses the value."""` | 0.0e+00 |
| test_phase1_adversarial.py::TestOracle::test_scalar_with_value_is_folded_by_engine_and_used_by_evaluator | `assert "m" not in str(eom[str(x.delta())][1])  # folded to the literal 2` | `assert "m" in str(eom[str(x.delta())][1])  # kept by name (D1)` | 0.0e+00 |
| test_s2_variation_fix.py::TestGyroscopicTerm::test_pinnedequation | `"(\\frac{d}{dt}(\\omega_{q})(-1)mll+Cross(\\omega_{q},c)(2)+Cross(q,e3)mgl(-1)+f)"` | `"(\\frac{d}{dt}(\\omega_{q})mll(-1)+Cross(\\omega_{q},c)(2)+Cross(q,e3)mgl(-1)+f)"` | 0.0e+00 |
| test_s2_variation_fix.py::TestAnisotropicInertiaOmegaForm::test_pinnedequation | `"(A\\frac{d}{dt}(\\omega_{q})(-1.0)+Cross(\\omega_{q},A\\omega_{q})(2.0)"` | `"(A\\frac{d}{dt}(\\omega_{q})(-1)+Cross(\\omega_{q},A\\omega_{q})(2)"` | 0.0e+00 |
| test_s2_variation_fix.py::TestRotatingSphere | `"(\\frac{d}{dt}(\\omega_{q})(-1)mLL""+Cross(q,Cross(e3,Cross(q,\\omega_{q})))\\OmegamLL""+Cross(q,Cross(\\omega_{q},Cross(e3,q)))(-1)\\OmegamLL""+Cross(Cross(Cross(\\omega_{q},q),e3),q)(-1)\\OmegamLL""+Cross(Cross(q,e3),Cross(\\omega_{q},q))(-1)\\OmegamLL""+Cross(\\omega_{q},Cross(q,Cross(e3,q)))\\OmegamLL(2.0)""+Cross(q,Cq)\\Omega\\OmegamLL+f)"` | `"(\\frac{d}{dt}(\\omega_{q})mLL(-1)+Cross(q,Cross(e3,Cross(q,\\omega_{q})))mLL\\Omega""+Cross(q,Cross(\\omega_{q},Cross(e3,q)))mLL\\Omega(-1)""+Cross(Cross(Cross(\\omega_{q},q),e3),q)mLL\\Omega(-1)""+Cross(Cross(q,e3),Cross(\\omega_{q},q))mLL\\Omega(-1)""+Cross(\\omega_{q},Cross(q,Cross(e3,q)))mLL\\Omega(2)+Cross(q,Cq)mLL\\Omega\\Omega+f)"` | 8.0e-17 |
| test_s2_variation_fix.py::TestRotatingSphere | `"(Cross(q,Cross(\\frac{d}{dt}(\\omega_{q}),q))(-1)mLL""+Cross(q,Cross(e3,Cross(\\omega_{q},q)))\\Omega(-1)mLL""+Cross(Cross(\\omega_{q},q),Cross(q,e3))\\OmegamLL""+Cross(q,Cross(Cross(\\omega_{q},q),e3))\\OmegamLL""+Cross(\\omega_{q},Cross(q,Cross(e3,q)))\\OmegamLL(2.0)""+Cross(q,Cross(\\omega_{q},Cross(e3,q)))\\Omega(-1)mLL""+Cross(q,Cq)\\Omega\\OmegamLL+f)"` | `"(Cross(q,Cross(\\frac{d}{dt}(\\omega_{q}),q))mLL(-1)""+Cross(q,Cross(e3,Cross(\\omega_{q},q)))mLL\\Omega(-1)""+Cross(Cross(\\omega_{q},q),Cross(q,e3))mLL\\Omega""+Cross(q,Cross(Cross(\\omega_{q},q),e3))mLL\\Omega""+Cross(\\omega_{q},Cross(q,Cross(e3,q)))mLL\\Omega(2)""+Cross(q,Cross(\\omega_{q},Cross(e3,q)))mLL\\Omega(-1)+Cross(q,Cq)mLL\\Omega\\Omega""+f)"` | 8.0e-17 |
| test_s2_variation_fix.py::TestRotatingFrame::test_pinnedequation | `"(\\frac{d}{dt}(\\omega_{q})(-1)mll+Cross(\\omega_{q},W0)mll(2.0)""+Cross(q,W0)mllDot(W0,q)(-1.0)+Cross(q,e3)mgl(-1)+f)"` | `"(\\frac{d}{dt}(\\omega_{q})mll(-1)+Cross(\\omega_{q},W0)mll(2)""+Cross(q,W0)mllDot(W0,q)(-1)+Cross(q,e3)mgl(-1)+f)"` | 0.0e+00 |
| test_s2_variation_fix.py::TestR3xS2VelocityCoupling::test_pinned_equations | `"(\\frac{d}{dt}(\\frac{d}{dt}(x))(-1)m+\\frac{d}{dt}(\\omega_{q})(-1)+F)"` | `"(\\frac{d}{dt}(\\frac{d}{dt}(x))m(-1)+\\frac{d}{dt}(\\omega_{q})(-1)+F)"` | 0.0e+00 |
| test_s2_variation_fix.py::TestR3xS2VelocityCoupling::test_pinned_equations | `"(\\frac{d}{dt}(\\omega_{q})(-1)mll+\\frac{d}{dt}(\\frac{d}{dt}(x))(-1)"` | `"(\\frac{d}{dt}(\\omega_{q})mll(-1)+\\frac{d}{dt}(\\frac{d}{dt}(x))(-1)"` | 0.0e+00 |
| test_s2_variation_fix.py::TestS2xSO3Couplings::test_omega_dot_Omega | `"(\\frac{d}{dt}(\\omega_{q})(-1)m+\\frac{d}{dt}(\\Omega_{R})(-1)k"` | `"(\\frac{d}{dt}(\\omega_{q})m(-1)+\\frac{d}{dt}(\\Omega_{R})k(-1)"` | 0.0e+00 |
| test_s2_variation_fix.py::TestS2xSO3Couplings::test_omega_dot_Omega | `"(J\\frac{d}{dt}(\\Omega_{R})(-1.0)+Cross(\\Omega_{R},J\\Omega_{R})(-1)""+\\frac{d}{dt}(\\omega_{q})(-1)k+Cross(\\Omega_{R},\\omega_{q})(-1)k+M)"` | `"(J\\frac{d}{dt}(\\Omega_{R})(-1)+Cross(\\Omega_{R},J\\Omega_{R})(-1)""+\\frac{d}{dt}(\\omega_{q})k(-1)+Cross(\\Omega_{R},\\omega_{q})k(-1)+M)"` | 0.0e+00 |
| test_s2_variation_fix.py::TestS2xSO3Couplings::test_omega_dot_body_axis | `"(\\frac{d}{dt}(\\omega_{q})(-1)m+Cross(\\Omega_{R},(R)'e3)k"` | `"(\\frac{d}{dt}(\\omega_{q})m(-1)+Cross(\\Omega_{R},(R)'e3)k"` | 0.0e+00 |
| test_s2_variation_fix.py::TestS2xSO3Couplings::test_omega_dot_body_axis | `"(J\\frac{d}{dt}(\\Omega_{R})(-1.0)+Cross(\\Omega_{R},J\\Omega_{R})(-1)""+Cross((R)'e3,\\omega_{q})(-1)k+M)"` | `"(J\\frac{d}{dt}(\\Omega_{R})(-1)+Cross(\\Omega_{R},J\\Omega_{R})(-1)""+Cross((R)'e3,\\omega_{q})k(-1)+M)"` | 0.0e+00 |
| test_scalar_adversarial.py::TestSpringMassExact::test_eom_exact | `assert str(eqn) == f"(m{DDX}(-1.0)+kx(-1.0)+f)"` | `assert str(eqn) == f"(m{DDX}(-1)+kx(-1)+f)"` | 0.0e+00 |
| test_scalar_adversarial.py::TestSpringMassExact::test_M_exact | `assert _strs(self.sf[self.key].M) == {DDX: "m(-1.0)"}` | `assert _strs(self.sf[self.key].M) == {DDX: "m(-1)"}` | 0.0e+00 |
| test_scalar_adversarial.py::TestSpringMassExact::test_f_exact | `assert str(self.sf[self.key].f) == "kx(-1.0)"` | `assert str(self.sf[self.key].f) == "kx(-1)"` | 0.0e+00 |
| test_scalar_adversarial.py::TestSpringMassExact::test_sign_convention_matches_s2_pendulum | `assert str(s2[str(xi)].M[DW]) == "I(-1)mll"assert str(self.sf[self.key].M[DDX]) == "m(-1.0)"` | `assert str(s2[str(xi)].M[DW]) == "Imll(-1)"assert str(self.sf[self.key].M[DDX]) == "m(-1)"` | 0.0e+00 |
| test_scalar_adversarial.py::TestSpringMassDamperExact::test_eom_exact | `assert str(eqn) == f"(m{DDX}(-1.0)+kx(-1.0)+f+c\\frac{{d}}{{dt}}(x)(-1))"` | `assert str(eqn) == f"(m{DDX}(-1)+kx(-1)+f+c\\frac{{d}}{{dt}}(x)(-1))"` | 0.0e+00 |
| test_scalar_adversarial.py::TestSpringMassDamperExact::test_M_exact | `assert _strs(self.sf[self.key].M) == {DDX: "m(-1.0)"}` | `assert _strs(self.sf[self.key].M) == {DDX: "m(-1)"}` | 0.0e+00 |
| test_scalar_adversarial.py::TestSpringMassDamperExact::test_f_exact | `assert str(self.sf[self.key].f) == "(kx(-1.0)+c\\frac{d}{dt}(x)(-1))"` | `assert str(self.sf[self.key].f) == "(kx(-1)+c\\frac{d}{dt}(x)(-1))"` | 0.0e+00 |
| test_scalar_adversarial.py::TestSpringMassDamperExact | `@pytest.mark.xfail(strict=True,reason="BUG: numeric coefficients mix int and float ((-1) and (-1.0)) in one equation",)` | `` | 0.0e+00 |
| test_scalar_adversarial.py::TestSpringMassDamperExact::test_numeric_coefficients_have_one_type | `''(-1.0)''. The S² pipeline strings use ''(-1)'' throughout."""` | `''(-1.0)''. The S² pipeline strings use ''(-1)'' throughout.Fixed in phase 5a (a float with an integer value folds to an int)."""` | 0.0e+00 |
| test_scalar_adversarial.py::TestDoubleSpringMassExact::test_x1_eom_exact | `assert str(eqn) == f"(m1{ddx1}(-1.0)+k1x1(-1.0)+k2x2+k2x1(-1.0)+f1)"` | `assert str(eqn) == f"(m1{ddx1}(-1)+k1x1(-1)+k2x2+k2x1(-1)+f1)"` | 0.0e+00 |
| test_scalar_adversarial.py::TestDoubleSpringMassExact::test_x2_eom_exact | `assert str(eqn) == f"(m2{ddx2}(-1.0)+k2x2(-1.0)+k2x1+f2)"` | `assert str(eqn) == f"(m2{ddx2}(-1)+k2x2(-1)+k2x1+f2)"` | 0.0e+00 |
| test_scalar_adversarial.py::TestDoubleSpringMassExact::test_x1_standard_form_exact | `assert _strs(eq.M) == {"\\frac{d}{dt}(\\frac{d}{dt}(x1))": "m1(-1.0)"}assert str(eq.f) == "(k1x1(-1.0)+k2x2+k2x1(-1.0))"` | `assert _strs(eq.M) == {"\\frac{d}{dt}(\\frac{d}{dt}(x1))": "m1(-1)"}assert str(eq.f) == "(k1x1(-1)+k2x2+k2x1(-1))"` | 0.0e+00 |
| test_scalar_adversarial.py::TestDoubleSpringMassExact::test_x2_standard_form_exact | `assert _strs(eq.M) == {"\\frac{d}{dt}(\\frac{d}{dt}(x2))": "m2(-1.0)"}assert str(eq.f) == "(k2x2(-1.0)+k2x1)"` | `assert _strs(eq.M) == {"\\frac{d}{dt}(\\frac{d}{dt}(x2))": "m2(-1)"}assert str(eq.f) == "(k2x2(-1)+k2x1)"` | 0.0e+00 |
| test_scalar_adversarial.py::TestCoupledScalarAccelerations::test_both_accelerations_in_M | `assert str(eom[str(Variation(x1))][1]) == f"(m{ddx1}(-1.0)+m{ddx2}(-1.0)+f1)"assert _strs(eq1.M) == {ddx1: "m(-1.0)", ddx2: "m(-1.0)"}assert _strs(eq2.M) == {ddx1: "m(-1.0)", ddx2: "m(-1.0)"}` | `assert str(eom[str(Variation(x1))][1]) == f"(m{ddx1}(-1)+m{ddx2}(-1)+f1)"assert _strs(eq1.M) == {ddx1: "m(-1)", ddx2: "m(-1)"}assert _strs(eq2.M) == {ddx1: "m(-1)", ddx2: "m(-1)"}` | 0.0e+00 |
| test_scalar_adversarial.py::TestVaryingInertia::test_mass_matrix_exact | `assert _strs(self.eq.M) == {DDX: "(m0(-1.0)+m1x(-1.0))"}` | `assert _strs(self.eq.M) == {DDX: "(m0(-1)+m1x(-1))"}` | 0.0e+00 |
| test_scalar_adversarial.py::TestVaryingInertia::test_f_is_mathematically_correct_but_uncombined | `"""The two ẋ² terms are present with coefficients +0.5 and −1.0."""assert str(self.eq.f) == ("(0.5m1\\frac{d}{dt}(x)\\frac{d}{dt}(x)+m1\\frac{d}{dt}(x)\\frac{d}{dt}(x)(-1.0))")` | `"""The two ẋ² terms are present with coefficients +0.5 and −1.0.The new printed form is the same expression simplified further (phase 5a mergesthe two terms into −0.5 m₁ ẋ²)."""assert str(self.eq.f) == ("m1\\frac{d}{dt}(x)\\frac{d}{dt}(x)(-0.5)")` | 0.0e+00 |
| test_scalar_adversarial.py::TestVaryingInertia | `@pytest.mark.xfail(strict=True,reason="BUG: like terms with a numeric factor buried inside the Mul chain do not combine",)` | `` | 0.0e+00 |
| test_scalar_adversarial.py::TestVaryingInertia::test_f_combines_like_terms | `"""½ m₁ ẋ² − m₁ ẋ² must collapse to a single term −½ m₁ ẋ²."""` | `"""½ m₁ ẋ² − m₁ ẋ² must collapse to a single term −½ m₁ ẋ².Fixed in phase 5a (like terms are compared on flattened Mul chains)."""` | 0.0e+00 |
| test_scalar_adversarial.py::TestVaryingInertia::test_x_squared_inertia_exact | `assert str(eqn) == f"(m{dx}x{dx}(-1.0)+mxx{DDX}(-1.0)+f)"` | `assert str(eqn) == f"(mx{dx}{dx}(-1)+mxx{DDX}(-1)+f)"` | 0.0e+00 |
| test_scalar_adversarial.py::TestVaryingInertia::test_x_squared_inertia_exact | `assert _strs(eq.M) == {DDX: "mxx(-1.0)"}` | `assert _strs(eq.M) == {DDX: "mxx(-1)"}` | 0.0e+00 |
| test_scalar_adversarial.py::TestElasticSphericalPendulum::test_l_eom_exact | `f"(m{DDL}(-1)+mlDot(\\omega_{{q}},\\omega_{{q}})+mgDot(q,e3)(-1)+kl(-1.0)+kl0+f_l)"` | `f"(m{DDL}(-1)+mlDot(\\omega_{{q}},\\omega_{{q}})+mgDot(q,e3)(-1)+kl(-1)+kl0+f_l)"` | 0.0e+00 |
| test_scalar_adversarial.py::TestElasticSphericalPendulum::test_l_standard_form_exact | `assert str(eq.f) == "(mlDot(\\omega_{q},\\omega_{q})+mgDot(q,e3)(-1)+kl(-1.0)+kl0)"` | `assert str(eq.f) == "(mlDot(\\omega_{q},\\omega_{q})+mgDot(q,e3)(-1)+kl(-1)+kl0)"` | 0.0e+00 |
| test_scalar_adversarial.py::TestElasticSphericalPendulum::test_xi_mass_and_input_exact | `assert _strs(eq.M) == {DW: "I(-1)mll"}` | `assert _strs(eq.M) == {DW: "Imll(-1)"}` | 0.0e+00 |
| test_scalar_adversarial.py::TestElasticSphericalPendulumQdotForm::test_l_eom_exact | `f"+m{DDL}(-1)+mgDot(q,e3)(-1)+kl(-1.0)+kl0+f_l)"` | `f"+m{DDL}(-1)+mgDot(q,e3)(-1)+kl(-1)+kl0+f_l)"` | 1.1e-16 |
| test_scalar_adversarial.py::TestElasticSphericalPendulumQdotForm | `@pytest.mark.xfail(strict=True,reason="F6 (deferred): the non-IBP pair +m l l̇ ω − m l l̇ ω is not combined, so f"" has three addends instead of two. The IBP Coriolis term −2 m l l̇ ω is present"" since F1.",)` | `` | 1.1e-16 |
| test_scalar_adversarial.py::TestElasticSphericalPendulumQdotForm::test_xi_f_exact_structure | `"""f must be −2 m l l̇ ω − m g l q×e₃ (two addends)."""` | `"""f must be −2 m l l̇ ω − m g l q×e₃ (two addends).Fixed in phase 5a (the pair +m l l̇ ω − m l l̇ ω now cancels)."""` | 1.1e-16 |
| test_scalar_adversarial.py::TestIntegrateByPartsScalarFactors::test_scalar_target_inside_dot_after_pull | `"(\\delta{l}(-1)Dot(\\frac{d}{dt}(q),w)+\\delta{l}(-1)Dot(q,\\frac{d}{dt}(w)))"` | `"(\\delta{l}Dot(\\frac{d}{dt}(q),w)(-1)+\\delta{l}Dot(q,\\frac{d}{dt}(w))(-1))"` | 0.0e+00 |
| test_scalar_adversarial.py::TestIntegrateByPartsScalarFactors::test_scalar_target_with_nonconstant_prefactor | `"(\\delta{l}(-1)\\frac{d}{dt}(a)Dot(q,w)""+\\delta{l}(-1)aDot(\\frac{d}{dt}(q),w)""+\\delta{l}(-1)aDot(q,\\frac{d}{dt}(w)))"` | `"(\\delta{l}\\frac{d}{dt}(a)Dot(q,w)(-1)+\\delta{l}aDot(\\frac{d}{dt}(q),w)(-1)""+\\delta{l}aDot(q,\\frac{d}{dt}(w))(-1))"` | 1.1e-16 |
| test_scalar_adversarial.py::TestNameCollisions::test_x_and_x1_coupled_exact | `assert str(eom[str(Variation(x))][1]) == f"(m{DDX}(-1.0)+kx1(-1)+f)"assert str(eom[str(Variation(x1))][1]) == f"(m{ddx1}(-1.0)+kx(-1)+f1)"assert _strs(sf[str(Variation(x))].M) == {DDX: "m(-1.0)"}assert _strs(sf[str(Variation(x1))].M) == {ddx1: "m(-1.0)"}` | `assert str(eom[str(Variation(x))][1]) == f"(m{DDX}(-1)+kx1(-1)+f)"assert str(eom[str(Variation(x1))][1]) == f"(m{ddx1}(-1)+kx(-1)+f1)"assert _strs(sf[str(Variation(x))].M) == {DDX: "m(-1)"}assert _strs(sf[str(Variation(x1))].M) == {ddx1: "m(-1)"}` | 0.0e+00 |
| test_scalar_adversarial.py::TestStandardFormHygiene::test_l_equation_mass_and_input | `assert str(self.eq_l.M[DDL]) == "mDot(Rrho,Rrho)(-1.0)"` | `assert str(self.eq_l.M[DDL]) == "mDot(Rrho,Rrho)(-1)"` | 1.4e-16 |
| test_scalar_adversarial.py::TestFixF1PipelineExact::test_elastic_pendulum_xi_equation_exact | `assert str(eqn) == (f"({DW}(-1)mll+\\omega_{{q}}ml{DL}(-2)+Cross(q,e3)mgl(-1)+f_q)")` | `assert str(eqn) == (f"({DW}mll(-1)+\\omega_{{q}}ml{DL}(-2)+Cross(q,e3)mgl(-1)+f_q)")` | 0.0e+00 |
| test_scalar_adversarial.py::TestFixF1PipelineExact::test_prescribed_length_s2_pendulum_exact | `assert str(eqn) == f"({DW}(-1)mll+\\omega_{{q}}ml{DL}(-2)+Cross(q,e3)mgl(-1)+f)"` | `assert str(eqn) == f"({DW}mll(-1)+\\omega_{{q}}ml{DL}(-2)+Cross(q,e3)mgl(-1)+f)"` | 0.0e+00 |
| test_scalar_adversarial.py::TestFixF1PipelineExact::test_r3_x_r_cross_inertia_exact | `assert str(eqn) == f"({DDX}(-1)m+e1m{DDL}(-1)+F)"assert _strs(sf[kx].M) == {DDL: "e1m(-1)", DDX: "I(-1)m"}` | `assert str(eqn) == f"({DDX}m(-1)+e1m{DDL}(-1)+F)"assert _strs(sf[kx].M) == {DDL: "e1m(-1)", DDX: "Im(-1)"}` | 0.0e+00 |
| test_simplify.py::TestSimplifyRecursive::test_chain_mul_numerics | `"""Mul(Mul(2, 3), a) → Mul(6, a)  (bottom-up combines inner first)."""` | `"""Mul(Mul(2, 3), a) → Mul(6, a)  (bottom-up combines inner first).Since phase 5a the folded number is placed last: Mul(a, 6)."""` | n/a |
| test_vbl.py::TestPointMassR3::test_pinned_eom | `assert str(self.eom[str(self.dx)][1]) == "(\\frac{d}{dt}(\\frac{d}{dt}(x))(-1)m+F)"` | `assert str(self.eom[str(self.dx)][1]) == "(\\frac{d}{dt}(\\frac{d}{dt}(x))m(-1)+F)"` | 0.0e+00 |
| test_vbl.py::TestSpringMassR::test_pinned_eom | `expected = "(m\\frac{d}{dt}(\\frac{d}{dt}(x))(-1.0)+kx(-1.0)+f)"` | `expected = "(m\\frac{d}{dt}(\\frac{d}{dt}(x))(-1)+kx(-1)+f)"` | 0.0e+00 |
| test_vbl.py::TestSpringMassR::test_dynamic_row_exact | `assert str(self.lin.A[row][str(self.dx)]) == "k(-1.0)"` | `assert str(self.lin.A[row][str(self.dx)]) == "k(-1)"` | 0.0e+00 (same output in a builder chain) |
| test_vbl.py::TestFreeRigidBody::test_pinned_eom | `expected = "(J\\frac{d}{dt}(\\Omega_{R})(-1.0)+Cross(\\Omega_{R},J\\Omega_{R})(-1)+M)"` | `expected = "(J\\frac{d}{dt}(\\Omega_{R})(-1)+Cross(\\Omega_{R},J\\Omega_{R})(-1)+M)"` | 0.0e+00 |
| test_vbl.py::TestRigidPendulumPivot::test_pinned_eom | `"(J\\frac{d}{dt}(\\Omega_{R})(-1.0)""+Cross(\\Omega_{R},J\\Omega_{R})(-1)""+Cross(\\rho,(R)'e3)mg(-1)""+M)"` | `"(J\\frac{d}{dt}(\\Omega_{R})(-1)+Cross(\\Omega_{R},J\\Omega_{R})(-1)""+Cross(\\rho,(R)'e3)mg(-1)+M)"` | 0.0e+00 |
| test_vbl.py::TestRigidPendulumPivot::test_gravity_block_exact | `assert str(self.lin.A[row][str(self.eta)]) == "Hat(\\rho)Hat((R)'e3)(-1)mg"` | `assert str(self.lin.A[row][str(self.eta)]) == "Hat(\\rho)Hat((R)'e3)mg(-1)"` | 0.0e+00 (same output in a builder chain) |
| test_vbl.py::TestSphericalPendulumS2::test_pinned_eom | `expected = "(\\frac{d}{dt}(\\omega_{q})(-1)mll+Cross(q,e3)mgl(-1)+f)"` | `expected = "(\\frac{d}{dt}(\\omega_{q})mll(-1)+Cross(q,e3)mgl(-1)+f)"` | 0.0e+00 |
| test_vbl.py::TestQuadrotor::test_pinned_eom | `"(\\frac{d}{dt}(\\frac{d}{dt}(x))(-1)m+e3mg(-1)+Re3f)"` | `"(\\frac{d}{dt}(\\frac{d}{dt}(x))m(-1)+e3mg(-1)+Re3f)"` | 0.0e+00 |
| test_vbl.py::TestQuadrotor::test_translational_row_exact | `assert str(self.lin.A[row][str(self.eta)]) == "RHat(e3)(-1)f"` | `assert str(self.lin.A[row][str(self.eta)]) == "RHat(e3)f(-1)"` | 0.0e+00 (same output in a builder chain) |
| test_vbl_adversarial.py::TestExactPins::test_scalar_damping_goes_to_velocity_block | `assert str(lin.A[row][str(x.delta())]) == "k(-1.0)"` | `assert str(lin.A[row][str(x.delta())]) == "k(-1)"` | 0.0e+00 |
| test_vbl_adversarial.py::TestErrorPathsHeld | `ids=["transpose_madd", "transpose_outer", "vee_madd", "vee_hat_vadd"],` | `ids=["transpose_outer", "vee_madd"],` | 0.0e+00 |
| test_vbl_adversarial.py::TestErrorPathsHeld | `` | `@pytest.mark.parametrize("E",[MVMul(Transpose(Hat(_q) + Hat(_om)), _om),Vee(Hat(_q + _om)),],ids=["transpose_madd", "vee_hat_vadd"],)def test_formerly_unsupported_structures_linearize_correctly(self, E):"""Phase 5a: (Hat(a) + Hat(b))ᵀ simplifies to −(Hat(a) + Hat(b)) (a sum of skewmatrices is skew) and Vee distributes over the sum, so these two used to raise andnow linearize. Not dropped: the blocks match finite differences."""eom, variables = {str(_xi): (_xi, E)}, SystemVariables(vectors=[_q])lin = linearize(eom, variables)assert finite_difference_check(lin, eom, variables, []) < 1e-6` | 0.0e+00 |
