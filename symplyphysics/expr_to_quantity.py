from sympy import Expr
from sympy.physics.units import Quantity
from sympy.physics.units.systems.si import SI

def convert(expr: Expr, quantity_name: str) -> Quantity:
  quantity_scale = SI._collect_factor_and_dimension(expr)
  result = Quantity(quantity_name)
  dimsys_SI = SI.get_dimension_system()
  dimsys_SI.set_quantity_dimension(result, quantity_scale[1])
  dimsys_SI.set_quantity_scale_factor(result, quantity_scale[0])
  return result
