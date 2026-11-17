# Models

Each page animates a system, draws it and names its variables, and shows the
kinetic and potential energy, the Lagrangian and virtual work, the equations
`compute_eom` derives from them, the code, and the checks the equations pass.
Five pages of systems with $n$ bodies also derive the equations once for every $n$
([sums over n bodies](../indexed.md)) and then substitute $n = 5$. The math on
every page is generated from the library and tested against it
(`scripts/model_pages.py`, `tests/test_model_pages.py`).

```{toctree}
:caption: Particles and pendulums
:maxdepth: 1

point_mass
spring_mass
planar_pendulum
double_planar_pendulum
cart_planar_pendulum
spherical_pendulum
particle_on_rotating_sphere
double_spherical_pendulum
n_link_spherical_chain
elastic_spherical_pendulum
cart_spherical_pendulum
```

```{toctree}
:caption: Rigid bodies
:maxdepth: 1

free_rigid_body
heavy_top
rigid_pendulum
double_rigid_pendulum
n_link_rigid_chain
two_rigid_bodies_spherical_joint
rigid_body
biped_3d
```

```{toctree}
:caption: Spacecraft
:maxdepth: 1

spacecraft_reaction_wheels
spacecraft_reaction_wheels_s1
spacecraft_solar_arrays
```

```{toctree}
:caption: Aerial vehicles
:maxdepth: 1

quadrotor
quadrotor_offset_com
quadrotor_suspended_load
quadrotor_elastic_load
quadrotor_flexible_cable
multi_quadrotor_point_mass
multi_quadrotor_rigid_load
multi_quadrotor_hose
```
