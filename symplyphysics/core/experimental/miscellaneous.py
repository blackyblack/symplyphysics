from typing import TypeVar, Any, Optional, Callable

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


def evaluate_or_global_fallback(evaluate: Optional[bool]) -> bool:
    """Fallbacks to the global value of `evaluate` if it is `None`, else returns the input."""

    if evaluate is not None:
        return evaluate

    return global_parameters.evaluate


def const(value: _T) -> Callable[..., _T]:
    """Returns a function that always returns the given `value`."""

    def closure(*_: Any) -> _T:
        return value

    return closure


__all__ = [
    "cacheit",
    "sympify_expr",
    "evaluate_or_global_fallback",
    "const",
]
