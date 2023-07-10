from sympy import Expr, sympify
from sympy.physics.units import Quantity as SymQuantity
from ..core.quantity_decorator import assert_equivalent_dimension
from ..core.symbols.quantities import Quantity


def convert_to(value: Quantity, target_unit: SymQuantity) -> Expr:
    """
    Convert ``value`` to its scale factor with all of its units represented as ``target_unit``.
    """
    target_quantity = Quantity(target_unit)
    assert_equivalent_dimension(value, value.dimension.name, "convert_to",
        target_quantity.dimension)
    return sympify(value.scale_factor) * (1 / sympify(target_quantity.scale_factor))
