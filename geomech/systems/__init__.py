"""Verified systems library: known mechanical systems, ready to derive, verify, linearize,
generate code for, simulate and visualize.

``get(name, **kwargs)`` (or the builder ``geomech.systems.<name>(**kwargs)``) returns a
``System`` holding the Lagrangian ``L``, the virtual work ``dW``, the configuration
``variables`` and ``inputs`` that ``compute_eom``/``to_standard_form`` take, default
numeric ``params``, a consistent initial state ``x0``, the ``reference`` the model comes
from and, for most systems, ``scene(layout)``, a ``geomech.viz.Scene`` of its bodies.
``list_systems()`` lists them. Builders take keyword arguments only and derive nothing.

Every system passes ``check_eom`` (≤ 1e-5), and the systems whose equations Lee et al. give
in closed form are compared with them at random admissible states; the tests store each
system's equations in ``tests/golden/`` (see ``docs/systems.md``).

Example::

    from geomech import check_eom, to_standard_form
    from geomech.systems import get

    sys = get("double_spherical_pendulum")
    eom = sys.derive()
    check_eom(sys.L, sys.dW, sys.variables, eom, inputs=sys.inputs)
    sf = to_standard_form(eom, sys.variables, sys.inputs)
"""

from geomech.systems._base import (
    MAX_SIZE,
    Reference,
    System,
    SystemInfo,
    compared_text,
    get,
    list_systems,
    names,
    source_text,
    systems_table,
)
from geomech.systems.aerial import (
    multi_quadrotor_point_mass,
    multi_quadrotor_rigid_load,
    quadrotor,
    quadrotor_elastic_load,
    quadrotor_flexible_cable,
    quadrotor_offset_com,
    quadrotor_suspended_load,
)
from geomech.systems.multibody import biped_3d, spacecraft_reaction_wheels
from geomech.systems.rigid import (
    double_rigid_pendulum,
    free_rigid_body,
    heavy_top,
    n_link_rigid_chain,
    rigid_body,
    rigid_pendulum,
    two_rigid_bodies_spherical_joint,
)
from geomech.systems.spheres import (
    cart_spherical_pendulum,
    double_spherical_pendulum,
    elastic_spherical_pendulum,
    n_link_spherical_chain,
    particle_on_rotating_sphere,
    point_mass,
    spherical_pendulum,
    spring_mass,
)

__all__ = [
    "System",
    "Reference",
    "SystemInfo",
    "MAX_SIZE",
    "compared_text",
    "get",
    "list_systems",
    "names",
    "source_text",
    "systems_table",
    "point_mass",
    "spring_mass",
    "spherical_pendulum",
    "particle_on_rotating_sphere",
    "double_spherical_pendulum",
    "n_link_spherical_chain",
    "elastic_spherical_pendulum",
    "cart_spherical_pendulum",
    "free_rigid_body",
    "heavy_top",
    "rigid_pendulum",
    "double_rigid_pendulum",
    "n_link_rigid_chain",
    "two_rigid_bodies_spherical_joint",
    "rigid_body",
    "spacecraft_reaction_wheels",
    "biped_3d",
    "quadrotor",
    "quadrotor_offset_com",
    "quadrotor_suspended_load",
    "quadrotor_elastic_load",
    "quadrotor_flexible_cable",
    "multi_quadrotor_point_mass",
    "multi_quadrotor_rigid_load",
]
