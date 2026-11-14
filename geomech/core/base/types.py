from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class ExprType(Enum):
    SCALAR = 1
    VECTOR = 2
    MATRIX = 3


#: Attribute strings accepted by ``attr=[...]`` on every leaf. ``"TangentVector"`` is
#: appended by ``TS2``/``TSO3`` themselves and sets no flag.
VALID_ATTRS = frozenset(
    {
        "Constant",
        "Zero",
        "Ones",
        "Identity",
        "UnitNorm",
        "SymmetricMatrix",
        "SkewSymmetricMatrix",
        "OrthogonalMatrix",
        "Manifold",
        "TangentVector",
    }
)


#: ``attr`` string → the flags it sets.
_ATTR_FLAGS = {
    "Constant": ("is_constant",),
    "Zero": ("is_zero", "is_constant"),
    "Ones": ("is_ones", "is_constant"),
    "Identity": ("is_identity", "is_constant"),
    "UnitNorm": ("is_unit_norm",),
    "SymmetricMatrix": ("is_symmetric",),
    "SkewSymmetricMatrix": ("is_skew_symmetric",),
    "OrthogonalMatrix": ("is_orthogonal",),
    "Manifold": ("is_manifold",),
    "TangentVector": (),  # marker only; no rule reads it
}


@dataclass(frozen=True)
class ExprFlags:
    """The structural facts the rewrite rules key on, fixed when a leaf is built.

    Frozen: a leaf's flags never change after construction (``dataclasses.replace`` makes
    a modified copy).
    """

    is_constant: bool = False
    is_zero: bool = False
    is_ones: bool = False
    is_unit_norm: bool = False
    is_identity: bool = False
    is_symmetric: bool = False
    is_skew_symmetric: bool = False
    is_orthogonal: bool = False
    is_numeric: bool = False
    is_manifold: bool = False

    @classmethod
    def from_attr_list(cls, attrs: list[str] | None, name: str | None = None) -> ExprFlags:
        """Build the flags for a leaf from its ``attr`` strings.

        An attribute outside ``VALID_ATTRS`` raises ``ValueError`` naming the symbol
        (*name*) and the attribute instead of being ignored; a bare string instead of a
        list raises ``TypeError``. Attribute names are case-sensitive.
        """
        if isinstance(attrs, str):
            where = f" for symbol '{name}'" if name is not None else ""
            raise TypeError(
                f"attr must be a list of strings, got str {attrs!r}{where}; write attr=[{attrs!r}]"
            )
        if not attrs:
            return cls()
        values = {}
        for a in attrs:
            if not isinstance(a, str) or a not in _ATTR_FLAGS:
                where = f" for symbol '{name}'" if name is not None else ""
                raise ValueError(
                    f"attr: unknown attribute {a!r}{where}; valid attributes are"
                    f" {sorted(VALID_ATTRS)}"
                )
            for flag in _ATTR_FLAGS[a]:
                values[flag] = True
        return cls(**values)


@dataclass(frozen=True)
class ManifoldInfo:
    tangent_vector_name: str
    variation_vector_name: str
