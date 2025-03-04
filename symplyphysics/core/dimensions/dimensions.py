from __future__ import annotations

from typing import Any
from sympy import Expr, S
from sympy.physics import units
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import dimsys_SI

from ..errors import UnitsError
from .collect_quantity import collect_quantity_factor_and_dimension
from .miscellaneous import is_any_dimension, is_number


class AnyDimension(Dimension):  # type: ignore[misc]
    # pylint: disable-next=signature-differs
    def __new__(cls) -> AnyDimension:
        return super().__new__(cls, "any_dimension")  # type: ignore[no-any-return]

    def _eval_nseries(self, _x: Any, _n: Any, _logx: Any, _cdir: Any) -> Any:
        pass


any_dimension = AnyDimension()


def assert_equivalent_dimension(
    arg: SymQuantity | Any | Dimension,
    param_name: str,
    func_name: str,
    expected_unit: SymQuantity | Dimension,
) -> None:
    """
    Asserts if the dimension of the argument matches the provided unit.

    Args:
        arg: Number, quantity, expression made of numbers and quantities, or dimension.
        param_name: Name of the parameter of the calling function.
        func_name: Name of the calling function.
        expected_unit: Expression or dimension which `arg` is compared to.

    Raises:
        TypeError: If `arg` is a number, but `expected_unit` is not dimensionless.
        UnitsError: If the dimensions don't match otherwise, or when the scale factor of `arg` is not a number.
    """

    if not isinstance(expected_unit, Dimension):
        expected_scale_factor, expected_unit = collect_quantity_factor_and_dimension(expected_unit)

        if is_any_dimension(expected_scale_factor) or isinstance(expected_unit, AnyDimension):
            return

    # HACK: this allows to treat angle type as dimensionless
    expected_unit = expected_unit.subs("angle", S.One)

    if not isinstance(arg, Dimension):
        (scale_factor, arg) = collect_quantity_factor_and_dimension(arg)

        if not is_number(scale_factor):
            # NOTE: this should probably raise `ValueError` or `TypeError`
            raise UnitsError(f"Argument '{param_name}' to function '{func_name}' should "
                f"not contain free symbols: '{scale_factor}'")

        if is_any_dimension(scale_factor) or isinstance(arg, AnyDimension):
            return

    # HACK: this allows to treat angle type as dimensionless
    arg = arg.subs("angle", S.One)

    if dimsys_SI.is_dimensionless(arg) and not dimsys_SI.is_dimensionless(expected_unit):
        # NOTE: this should probably be `UnitsError`
        raise TypeError(f"Argument '{param_name}' to function '{func_name}'"
            f" is Number but '{expected_unit}' is not dimensionless")

    if not dimsys_SI.equivalent_dims(arg, expected_unit):
        raise UnitsError(f"Argument '{param_name}' to function '{func_name}' must "
            f"be in units equivalent to '{expected_unit.name}', got {arg.name}")


def print_dimension(dimension: Dimension) -> str:
    """Returns the prettified name of ``dimension``."""

    return "dimensionless" if dimsys_SI.is_dimensionless(dimension) else str(dimension.name)


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
    # re-exports
    "Dimension",
    "dimsys_SI",

    # locals
    "AnyDimension",
    "any_dimension",
    "assert_equivalent_dimension",
    "print_dimension",
    "dimension_to_si_unit",
]
