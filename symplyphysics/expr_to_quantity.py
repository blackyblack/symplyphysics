from sympy import Expr
from sympy.core.singleton import S
from sympy.functions import exp
from sympy.physics.units import Quantity
from sympy.physics.units.systems.si import SI

def expr_to_quantity(expr: Expr, quantity_name: str) -> Quantity:
    quantity_scale = SI._collect_factor_and_dimension(expr)
    dimension = quantity_scale[1]
    if isinstance(dimension, exp):
        dimension_args = dimension.nargs.args[0]
        # power of the exponent should be dimensionless
        assert dimension_args == S.One
        dimension = S.One
    result = Quantity(quantity_name)
    dimsys_SI = SI.get_dimension_system()
    dimsys_SI.set_quantity_dimension(result, dimension)
    dimsys_SI.set_quantity_scale_factor(result, quantity_scale[0])
    return result
