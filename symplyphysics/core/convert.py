from typing import Any
from sympy import Expr, S, sympify
from sympy.physics.units import Quantity as SymQuantity

from .dimensions import assert_equivalent_dimension, dimension_to_si_unit
from .symbols.quantities import Quantity


def convert_to(value: Expr, target_unit: Expr) -> Expr:
    """
    Convert ``value`` to its scale factor with ``value`` unit represented as ``target_unit``.
    """
    if not isinstance(value, SymQuantity):
        value = Quantity(value)
    if not isinstance(target_unit, SymQuantity):
        target_unit = Quantity(target_unit)

    assert_equivalent_dimension(value, value.dimension.name, "convert_to", target_unit.dimension)
    return sympify(value.scale_factor) * (1 / sympify(target_unit.scale_factor))


def convert_to_float(value: Expr) -> float:
    return float(convert_to(value, S.One))


def convert_to_si(value: Expr | float) -> Expr:
    if not isinstance(value, SymQuantity):
        value = Quantity(value)

    unit = dimension_to_si_unit(value.dimension)
    return convert_to(value, unit)


def evaluate_quantity(quantity: Expr, **kwargs: Any) -> Quantity:
    if not isinstance(quantity, SymQuantity):
        quantity = Quantity(quantity)

    scale_factor_ = quantity.scale_factor.evalf(**kwargs)
    dimension = quantity.dimension
    return Quantity(scale_factor_, dimension=dimension)


def evaluate_expression(expr: Expr, evaluate: bool = False, **kwargs: Any) -> Expr:
    for qty in expr.atoms(SymQuantity):
        si_value = convert_to_si(qty)
        if evaluate:
            si_value = si_value.evalf(**kwargs)
        expr = expr.subs(qty, si_value)

    return expr


__all__ = [
    "convert_to",
    "convert_to_float",
    "convert_to_si",
    "evaluate_quantity",
    "evaluate_expression",
]
