"""System variables for Lagrangian dynamics."""

from __future__ import annotations

from dataclasses import dataclass, field


def check_group_is_sequence(group_name: str, group: object) -> None:
    """A ``SystemVariables`` group must be a list (or tuple) of variables, not one leaf.

    A bare leaf (``vectors=x``) or a one-shot iterator would otherwise fail with Python's
    own ``'Vector' object is not iterable`` or be consumed by the first pass over it.
    """
    if not isinstance(group, (list, tuple)):
        hint = f"; write {group_name}=[{group}]" if hasattr(group, "type") else ""
        raise TypeError(
            f"SystemVariables.{group_name} must be a list of configuration variables, got"
            f" {type(group).__name__} {str(group)!r}{hint}"
        )


@dataclass
class SystemVariables:
    """Configuration variables that define the degrees of freedom.

    scalars:  scalar configuration variables (Scalar)
    vectors:  vector configuration variables (Vector, S2)
    matrices: matrix configuration variables (SO3)

    Each group is a list (or tuple); a bare leaf raises ``TypeError``. An empty
    ``SystemVariables()`` can be built (``check_eom`` accepts it), but ``compute_eom``
    raises ``ValueError`` for it: a system with no degree of freedom has no equation.
    """

    scalars: list = field(default_factory=list)
    vectors: list = field(default_factory=list)
    matrices: list = field(default_factory=list)

    def __post_init__(self):
        for group_name in ("scalars", "vectors", "matrices"):
            check_group_is_sequence(group_name, getattr(self, group_name))
