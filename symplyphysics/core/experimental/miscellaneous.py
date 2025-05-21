from typing import TypeVar, Any, Optional

from sympy import cacheit as sym_cacheit, Expr, sympify
from sympy.core.parameters import global_parameters

_T = TypeVar("_T")


def cacheit(f: _T) -> _T:
    """A typed version of `sympy.cacheit`."""

    return sym_cacheit(f)


def sympify_expr(value: Any) -> Expr:
    result = sympify(value, strict=True)

    if not isinstance(result, Expr):
        raise TypeError(f"Expected '{result}' to be Expr.")

    return result


def set_evaluate(evaluate: Optional[bool]) -> bool:
    if evaluate is not None:
        return evaluate

    return global_parameters.evaluate


__all__ = [
    "cacheit",
    "sympify_expr",
    "set_evaluate",
]
