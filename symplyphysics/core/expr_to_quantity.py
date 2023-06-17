from sympy import Expr, S, exp
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.vector import VectorAdd
from ..core.coordinate_systems.coordinate_systems import CoordinateSystem
from ..core.vectors.vectors import Vector, vector_from_sympy_vector
from ..core.symbols.quantities import Quantity, collect_factor_and_dimension


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


def expr_to_vector(expr: Expr, coordinate_system: CoordinateSystem = None) -> Vector:
    if isinstance(expr, Mul):
        expr = expr.expand()
    if isinstance(expr, Add):
        expr = VectorAdd(expr)
    vector = vector_from_sympy_vector(expr, coordinate_system)
    components = []
    for c in vector.components:
        components.append(expr_to_quantity(c))
    return Vector(components, coordinate_system)
