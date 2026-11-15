| Example | n | rank | max Re λ | min ζ | e₀ pos. (m) | e₀ att. (rad) | t_s | peak ‖u − u_eq‖∞ | drift | converged |
|---|---|---|---|---|---|---|---|---|---|---|
| `spherical_pendulum_inverted` | 4 | 4 | -3.17 | 1.00 | — | 0.5 | 1.39 s | 9.68 N·m | 7.0e-15 | yes |
| `double_spherical_pendulum_inverted` | 8 | 8 | -2.55 | 1.00 | — | 0.312 | 1.60 s | 6.33 N·m | 9.9e-15 | yes |
| `n_link_spherical_chain_inverted` | 12 | 12 | -2.17 | 1.00 | — | 0.305 | 1.88 s | 7.72 N·m | 1.5e-14 | yes |
| `elastic_spherical_pendulum_inverted` | 6 | 6 | -0.559 | 0.09 | 0.1 | 0.342 | 4.38 s | 0.527 N, 5.29 N·m | 7.1e-15 | yes |
| `cart_planar_pendulum_inverted` | 4 | 4 | -1.51 | 0.79 | 0.3 | 0.175 | 3.01 s | 21.4 N | 1.1e-16 | yes |
| `rigid_pendulum_inverted` | 6 | 6 | -3.16 | 1.00 | — | 0.436 | 1.22 s | 3.8 N·m | 4.1e-14 | yes |
| `double_rigid_pendulum_inverted` | 12 | 12 | -2.97 | 1.00 | — | 0.315 | 1.50 s | 4.32 N·m | 4.1e-14 | yes |
| `n_link_rigid_chain_inverted` | 18 | 18 | -2.42 | 1.00 | — | 0.306 | 1.57 s | 5.87 N·m | 5.1e-14 | yes |
| `free_rigid_body_attitude` | 6 | 6 | -1.66 | 0.78 | — | 1.92 | 2.27 s | 4.96 N·m | 1.1e-13 | yes |
| `spacecraft_reaction_wheels_slew` † | 9 | 6 | -0.371 | 0.72 | — | 1.05 | 10.15 s | 0.191 N·m | 7.0e-14 | yes |
| `spacecraft_reaction_wheels_s1_slew` † | 9 | 6 | -0.371 | 0.72 | — | 1.05 | 10.15 s | 0.191 N·m | 7.0e-14 | yes |
| `spacecraft_solar_arrays_slew` † | 13 | 10 | -0.281 | 0.73 | — | 1.44 | 13.52 s | 0.191 N·m | 1.3e-13 | yes |
| `biped_3d_balance` | 12 | 12 | -0.882 | 0.21 | — | 0.145 | 1.87 s | 20.4 N·m | 6.9e-14 | yes |
| `rigid_body_pose` | 12 | 12 | -0.924 | 0.73 | 1.22 | 1.57 | 3.93 s | 3.16 N, 4.97 N·m | 5.6e-14 | yes |
| `quadrotor_hover` | 12 | 12 | -0.924 | 0.71 | 1.87 | 1.05 | 3.78 s | 1.58 N, 4.76 N·m | 9.5e-14 | yes |
| `quadrotor_offset_com_hover` | 12 | 12 | -0.923 | 0.68 | 1.22 | 0.785 | 3.64 s | 1.53 N, 3.23 N·m | 1.3e-13 | yes |
| `quadrotor_suspended_load_hover` | 16 | 16 | -0.596 | 0.36 | 0.781 | 0.92 | 4.49 s | 0.761 N, 2.33 N·m | 7.3e-14 | yes |
| `quadrotor_elastic_load_hover` | 18 | 18 | -0.596 | 0.04 | 0.473 | 0.406 | 4.23 s | 2.51 N, 0.752 N·m | 9.7e-14 | yes |
| `quadrotor_flexible_cable_hover` | 24 | 24 | -0.166 | 0.01 | 0.374 | 0.203 | 21.70 s | 0.12 N, 0.808 N·m | 1.4e-13 | yes |
| `multi_quadrotor_point_mass_hover` | 36 | 36 | -0.45 | 0.48 | 0.412 | 0.139 | 6.79 s | 0.221 N, 0.276 N·m | 1.5e-13 | yes |
| `multi_quadrotor_rigid_load_hover` | 42 | 42 | -0.438 | 0.17 | 0.412 | 0.281 | 6.81 s | 0.143 N, 0.321 N·m | 1.3e-13 | yes |
| `two_rigid_bodies_spherical_joint_hold` | 18 | 18 | -1.03 | 0.74 | 0.412 | 0.554 | 3.64 s | 0.967 N, 1.11 N·m | 9.5e-14 | yes |

† Uncontrollable modes with eigenvalue 0 (conserved quantities); LQR on the controllable subspace.

**Not stabilizable** (uncontrollable modes listed; LQR on the controllable subspace):

| Example | n | rank | uncontrollable eigenvalues | max Re λ (controlled) | e₀ att. (rad) | t_s | converged |
|---|---|---|---|---|---|---|---|
| `cart_spherical_pendulum_inverted_in_plane` | 6 | 4 | -3.5, +3.5 | -1.51 | 0.174 | 3.01 s | yes (in-plane start only) |
| `cart_spherical_pendulum_out_of_plane` | 6 | 4 | -3.5, +3.5 | -1.51 | 0.174 | — | diverged |
