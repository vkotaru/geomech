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


@dataclass
class ExprFlags:
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
        flags = cls()
        for a in attrs:
            match a:
                case "Constant":
                    flags.is_constant = True
                case "Zero":
                    flags.is_zero = True
                    flags.is_constant = True
                case "Ones":
                    flags.is_ones = True
                    flags.is_constant = True
                case "Identity":
                    flags.is_identity = True
                    flags.is_constant = True
                case "UnitNorm":
                    flags.is_unit_norm = True
                case "SymmetricMatrix":
                    flags.is_symmetric = True
                case "SkewSymmetricMatrix":
                    flags.is_skew_symmetric = True
                case "OrthogonalMatrix":
                    flags.is_orthogonal = True
                case "Manifold":
                    flags.is_manifold = True
                case "TangentVector":
                    pass  # marker only; no rule reads it
                case _:
                    where = f" for symbol '{name}'" if name is not None else ""
                    raise ValueError(
                        f"attr: unknown attribute {a!r}{where}; valid attributes are"
                        f" {sorted(VALID_ATTRS)}"
                    )
        return flags


@dataclass
class ManifoldInfo:
    tangent_vector_name: str
    variation_vector_name: str
