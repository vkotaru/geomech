_OPERAND_MAX = 80


def _operand_str(expr):
    """``str(expr)`` quoted for an error message, shortened past ``_OPERAND_MAX`` characters."""
    s = str(expr)
    if len(s) > _OPERAND_MAX:
        s = s[: _OPERAND_MAX - 3] + "..."
    return "'" + s + "'"


def _short_repr(val):
    """``repr(val)`` for an error message, shortened past ``_OPERAND_MAX`` characters."""
    s = repr(val)
    return s if len(s) <= _OPERAND_MAX else s[: _OPERAND_MAX - 3] + "..."


def _type_name(t):
    return getattr(t, "name", str(t))


class ExpressionMismatchError(Exception):
    """Operand types do not fit the operation.

    Built from the operation name and the two operand types; when the operands themselves
    are passed (``lhs``/``rhs``) the message names them, so that the failing node can be
    found in a long Lagrangian:

        Dot cannot combine SCALAR 'm' and VECTOR 'x'

    ``ExpressionMismatchError.unary(op, expected, expr)`` builds the one-operand form
    (``Hat expects a VECTOR, got SCALAR 'm'``); ``msg`` replaces the message outright.
    """

    def __init__(self, operation="", ltype=None, rtype=None, lhs=None, rhs=None, *, msg=None):
        if msg is None:
            if ltype is not None and rtype is not None and (lhs is not None or rhs is not None):
                lpart = _type_name(ltype) + (" " + _operand_str(lhs) if lhs is not None else "")
                rpart = _type_name(rtype) + (" " + _operand_str(rhs) if rhs is not None else "")
                msg = f"{operation} cannot combine {lpart} and {rpart}"
            elif ltype is not None and rtype is not None:
                msg = (
                    operation
                    + " cannot be performed between "
                    + _type_name(ltype)
                    + " and "
                    + _type_name(rtype)
                )
            else:
                msg = (operation + ": " if operation else "") + "operation cannot be performed"
        super().__init__(msg)

    @classmethod
    def unary(cls, operation, expected, expr):
        """``{operation} expects a {expected}, got {type} '{expr}'`` for a one-operand node."""
        return cls(
            msg=(
                f"{operation} expects a {_type_name(expected)},"
                f" got {_type_name(expr.type)} {_operand_str(expr)}"
            )
        )


class SizeMismatchError(Exception):
    """Operand sizes do not match for the operation."""

    def __init__(self, operation="", expected=None, got=None):
        if expected is not None and got is not None:
            msg = operation + " size mismatch: expected " + str(expected) + ", got " + str(got)
        else:
            msg = operation + " operand sizes do not match"
        super().__init__(msg)


class AlgebraicError(Exception):
    """Expression cannot be algebraically factored as requested.

    Raised when an operation is mathematically impossible, not merely
    unimplemented — e.g. extracting a linear vector coefficient from
    an outer product (rank-1 matrix) that does not factor as C * q.
    """

    def __init__(self, msg=""):
        super().__init__(msg)


class NonLinearError(Exception):
    """Query appears nonlinearly (degree > 1) in the expression.

    Raised when linear coefficient extraction encounters the query
    in both children of a binary node — e.g. Mul(f(q), g(q)).
    """

    def __init__(self, msg=""):
        super().__init__(msg)
