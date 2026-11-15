geomech
=======

.. automodule:: geomech

``import geomech`` re-exports the public names of the subpackages below, so
``from geomech import compute_eom, Dot, S2`` works; each name is documented on the page of
the module that defines it.

==================================================  ===========================================
Page                                                Contents
==================================================  ===========================================
:doc:`core-base`                                    types, flags, leaves ``Scalar``/``Vector``/``Matrix``, ``S1``, ``S2``, ``SO3``
:doc:`core-operations`                              the node classes (sums, products, geometry, calculus)
:doc:`core-transformations`                         ``expand``, ``pull``, ``simplify``, manifold rules
:doc:`core-math`                                    ``extract_linear_coeff``, ``collect``, integration by parts
:doc:`dynamics`                                     ``compute_eom``, ``to_standard_form``, ``linearize``
:doc:`numeric`                                      the numeric oracle (``evaluate``, ``check_eom``)
:doc:`codegen`                                      ``StateLayout``, ``lambdify``, ``generate_module``
:doc:`simulation`                                   ``simulate``, ``accelerations``, ``lqr``
:doc:`viz`                                          ``Scene``, ``Report``
:doc:`systems`                                      the verified systems library (``get``, ``list_systems``)
:doc:`utils`                                        printing, tracing, invariants, error types
==================================================  ===========================================
