from typing import Any
from sympy import Expr, S
from sympy.physics.units import Dimension


def is_any_dimension(factor: Expr) -> bool:
    """
    Checks if ``factor`` is `0`, `±Inf`, or `NaN`, which can have any dimension due to their
    absorbing nature.
    """

    return factor in (S.Zero, S.Infinity, S.NegativeInfinity, S.NaN)


def is_number(value: Any) -> bool:
    """Checks if ``value`` is a (complex) number."""

    try:
        complex(value)
    except (TypeError, ValueError):
        return False

    return True


dimensionless = Dimension(S.One)
"""Alias for `Dimension(1)`."""

__all__ = [
    "is_any_dimension",
    "is_number",
    "dimensionless",
]
