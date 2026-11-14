geomech.systems
===============

.. automodule:: geomech.systems
   :ignore-module-all:
   :no-members:

.. autoclass:: geomech.systems.System
   :members:

.. autoclass:: geomech.systems.Reference

.. autoclass:: geomech.systems.SystemInfo

.. autofunction:: geomech.systems.get

.. autofunction:: geomech.systems.list_systems

.. autofunction:: geomech.systems.names

.. autofunction:: geomech.systems.systems_table

.. autofunction:: geomech.systems.source_text

Points, springs and pendulums: ``geomech.systems.spheres``
-----------------------------------------------------------

.. automodule:: geomech.systems.spheres
   :members: point_mass, spring_mass, spherical_pendulum, particle_on_rotating_sphere,
             double_spherical_pendulum, n_link_spherical_chain, elastic_spherical_pendulum,
             cart_spherical_pendulum

Rigid bodies: ``geomech.systems.rigid``
---------------------------------------

.. automodule:: geomech.systems.rigid
   :members: free_rigid_body, heavy_top, rigid_pendulum, double_rigid_pendulum,
             two_rigid_bodies_spherical_joint, rigid_body

Aerial vehicles: ``geomech.systems.aerial``
-------------------------------------------

.. automodule:: geomech.systems.aerial
   :members: quadrotor, quadrotor_offset_com, quadrotor_suspended_load,
             quadrotor_elastic_load, quadrotor_flexible_cable, multi_quadrotor_point_mass,
             multi_quadrotor_rigid_load, add_quadrotor
