from functools import reduce
from typing import Optional, Iterable
from sympy import Expr, S, exp, Add, Mul, Pow, sympify
from sympy.physics.units.systems.si import SI
from sympy.physics.units.util import _get_conversion_matrix_for_expr
from sympy.vector import VectorAdd, VectorMul, Vector as SymVector
from sympy.physics.units import Quantity as SymQuantity
from ..core.coordinate_systems.coordinate_systems import CoordinateSystem
from ..core.vectors.vectors import Vector, vector_from_sympy_vector
from ..core.symbols.quantities import Quantity, collect_factor_and_dimension


# HACK: copy of SymPy convert_to() with fixed prefixed units processing. It
# also removes units from the converted result, so one does not need
# to subs() all units to get real value.
def convert_to(expr: Expr, target_units: Iterable[SymQuantity] | SymQuantity) -> Expr:
    target_units_list: list[SymQuantity]
    if not isinstance(target_units, Iterable):
        target_units_list = [target_units]
    else:
        target_units_list = list(target_units)

    if isinstance(expr, Add):
        return Add.fromiter(convert_to(i, target_units_list) for i in expr.args)

    expr = sympify(expr)
    target_units_list = sympify(target_units_list)

    if not isinstance(expr, SymQuantity) and expr.has(SymQuantity):
        expr = expr.replace(lambda x: isinstance(x, SymQuantity),
            lambda x: convert_to(x, target_units_list))

    def get_total_scale_factor(expr: Expr) -> Expr:
        if isinstance(expr, Mul):
            return reduce(lambda x, y: x * y, [get_total_scale_factor(i) for i in expr.args])
        if isinstance(expr, Pow):
            return get_total_scale_factor(expr.base)**expr.exp
        if isinstance(expr, SymQuantity):
            return SI.get_quantity_scale_factor(expr)
        return expr

    depmat = _get_conversion_matrix_for_expr(expr, target_units_list, SI)
    if depmat is None:
        return expr

    expr_scale_factor = get_total_scale_factor(expr)
    vals: list[Expr] = []
    for u, p in zip(target_units_list, depmat):
        s = get_total_scale_factor(u)
        vals.append((1 / s)**p)
    return expr_scale_factor * Mul.fromiter(vals)


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
