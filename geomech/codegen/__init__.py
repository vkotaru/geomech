"""Code generation: numeric functions from the symbolic output of the pipeline.

``StateLayout`` maps the symbols of a system to flat arrays; ``to_source``/``lambdify``
turn one expression into Python source or a callable; ``generate_module`` writes a module
with the blocks of a ``compute_eom``, ``to_standard_form`` or ``linearize`` result;
``energy_function`` generates the total energy of a Lagrangian. Two targets: ``"numpy"``
(default, no dependency beyond numpy) and ``"casadi"`` (optional extra
``geomech[casadi]``; imported only when that target is used). See ``docs/codegen.md``.
"""

from geomech.codegen.energy import energy_expr, energy_function, energy_source, momenta
from geomech.codegen.layout import Slot, StateLayout, collect_params
from geomech.codegen.module import generate_module, load_module, module_source
from geomech.codegen.source import lambdify, to_source

__all__ = [
    "StateLayout",
    "Slot",
    "collect_params",
    "to_source",
    "lambdify",
    "generate_module",
    "module_source",
    "load_module",
    "energy_expr",
    "energy_function",
    "energy_source",
    "momenta",
]
