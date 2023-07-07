from sympy import Expr, sympify
from sympy.physics.units import Quantity as SymQuantity
from ..core.quantity_decorator import assert_equivalent_dimension
from ..core.symbols.quantities import collect_factor_and_dimension


def convert_to(value: SymQuantity, target_unit: SymQuantity) -> Expr:
    """
    Convert ``value`` to its scale factor with all of its units represented as ``target_unit``.
    """
    (target_scale, target_dim) = collect_factor_and_dimension(target_unit)
    (value_scale, _) = collect_factor_and_dimension(value)
    assert_equivalent_dimension(value, value.dimension.name, "convert_to", target_dim)
    return sympify(value_scale) * (1 / sympify(target_scale))
