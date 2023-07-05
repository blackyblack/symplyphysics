from functools import reduce
from typing import Optional
from sympy import Expr, S, exp, Add, Mul, Pow
from sympy.physics.units import Quantity as SymQuantity, Dimension
from sympy.physics.units.systems.si import SI
from sympy.vector import VectorAdd, VectorMul, Vector as SymVector

from ..core.quantity_decorator import assert_equivalent_dimension
from ..core.coordinate_systems.coordinate_systems import CoordinateSystem
from ..core.vectors.vectors import Vector, vector_from_sympy_vector
from ..core.symbols.quantities import Quantity, collect_factor_and_dimension


def _get_total_scale_factor(expr: Expr) -> Expr:
    if isinstance(expr, Mul):
        return reduce(lambda x, y: x * y, [_get_total_scale_factor(i) for i in expr.args])
    if isinstance(expr, Pow):
        return _get_total_scale_factor(expr.base)**expr.exp
    if isinstance(expr, SymQuantity):
        return SI.get_quantity_scale_factor(expr)
    return expr


def convert_to(value: SymQuantity, target_unit: SymQuantity) -> SymQuantity:
    """
    Convert ``value`` to the same Quantity with all of its units represented as ``target_unit``.
    """
    target_dim = Dimension(SI.get_dimensional_expr(target_unit))
    assert_equivalent_dimension(value, value.dimension.name, "convert_to", target_dim)
    expr_scale_factor = _get_total_scale_factor(value)
    s = _get_total_scale_factor(target_unit)
    return expr_scale_factor * (1 / s)


def expr_to_quantity(expr: Expr) -> Quantity:
    quantity_scale = collect_factor_and_dimension(expr)
    dimension = quantity_scale[1]
    # HACK: this allows to treat angle type as dimensionless
    dimension = dimension.subs("angle", S.One)
    if isinstance(dimension, exp):
        dimension_args = dimension.nargs.args[0]
        # power of the exponent should be dimensionless
        assert dimension_args == S.One
        dimension = S.One
    return Quantity(quantity_scale[0], dimension=dimension)


def expr_to_vector(expr: Expr, coordinate_system: Optional[CoordinateSystem] = None) -> Vector:
    if isinstance(expr, Mul):
        expr = VectorMul(*expr.args)
    if isinstance(expr, Add):
        expr = VectorAdd(*expr.args)
    if not isinstance(expr, SymVector):
        raise TypeError(f"Expression cannot be converted to SymPy Vector: {str(expr)}")
    vector = vector_from_sympy_vector(expr, coordinate_system)
    components: list[Quantity] = []
    for c in vector.components:
        components.append(expr_to_quantity(c))
    return Vector(components, coordinate_system)
