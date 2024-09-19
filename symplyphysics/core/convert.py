from sympy import Expr, sympify, S
from sympy.physics import units
from sympy.physics.units.systems.si import dimsys_SI
from sympy.physics.units.definitions import dimension_definitions

from .dimensions import assert_equivalent_dimension
from .symbols.quantities import Quantity


def convert_to(value: Expr, target_unit: Expr) -> Expr:
    """
    Convert ``value`` to its scale factor with ``value`` unit represented as ``target_unit``.
    """
    value_quantity = value if isinstance(value, Quantity) else Quantity(value)
    target_quantity = target_unit if isinstance(target_unit, Quantity) else Quantity(target_unit)
    assert_equivalent_dimension(value_quantity, value_quantity.dimension.name, "convert_to",
        target_quantity.dimension)
    return sympify(value_quantity.scale_factor) * (1 / sympify(target_quantity.scale_factor))


def convert_to_float(value: Expr) -> float:
    return float(convert_to(value, S.One))


_si_conversions: dict[units.Dimension, Expr] = {
    dimension_definitions.angle: S.One,
    dimension_definitions.length: units.meter,
    dimension_definitions.mass: units.kilogram,
    dimension_definitions.time: units.second,
    dimension_definitions.current: units.ampere,
    dimension_definitions.temperature: units.kelvin,
    dimension_definitions.amount_of_substance: units.mole,
    dimension_definitions.luminous_intensity: units.candela,
}


def convert_to_si(value: Expr | float) -> Expr:
    quantity = value if isinstance(value, Quantity) else Quantity(value)
    dependencies = dimsys_SI.get_dimensional_dependencies(quantity.dimension)
    unit = S.One
    for dimension, power in dependencies.items():
        unit *= _si_conversions[dimension]**power
    return convert_to(quantity, unit)


__all__ = [
    "convert_to",
    "convert_to_float",
    "convert_to_si",
]
