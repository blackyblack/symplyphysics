from typing import Any

from sympy import Expr, S
from sympy.physics import units
from sympy.physics.units import Dimension
from sympy.physics.units.systems.si import dimsys_SI


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

_si_conversions = {
    units.length: units.meter,
    units.mass: units.kilogram,
    units.time: units.second,
    units.current: units.ampere,
    units.temperature: units.kelvin,
    units.amount_of_substance: units.mole,
    units.luminous_intensity: units.candela,
}


def dimension_to_si_unit(dimension: Dimension) -> Expr:
    """Converts ``dimension`` to the corresponding SI unit."""

    si_unit = S.One

    dependencies = dimsys_SI.get_dimensional_dependencies(dimension)
    for dim, n in dependencies.items():
        si_unit *= _si_conversions.get(dim, S.One)**n

    return si_unit


__all__ = [
    "is_any_dimension",
    "is_number",
    "dimensionless",
    "dimension_to_si_unit",
]
