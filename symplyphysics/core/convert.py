from sympy import Expr, sympify, S
# from sympy.physics.units import Quantity as SymQuantity

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
