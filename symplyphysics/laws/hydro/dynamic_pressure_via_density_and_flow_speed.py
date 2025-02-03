"""
Dynamic pressure from speed
===========================

**Dynamic pressure**, sometimes called **velocity pressure**, is a physical quantity denoting the
pressure caused by a flowing fluid.

**Notes:**

#. Many authors define this quantity only for *incompressible flows*, but others extend it for
   compressible flows as well.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Dynamic_pressure>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

dynamic_pressure = symbols.dynamic_pressure
"""
:symbols:`dynamic_pressure` of the fluid. 
"""

density = symbols.density
"""
:symbols:`density` of the fluid.
"""

flow_speed = symbols.flow_speed
"""
:symbols:`flow_speed` of the fluid.
"""

law = Eq(dynamic_pressure, density * flow_speed**2 / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(density_=density, velocity_=flow_speed)
@validate_output(dynamic_pressure)
def calculate_pressure(density_: Quantity, velocity_: Quantity) -> Quantity:
    result_pressure_expr = solve(law, dynamic_pressure, dict=True)[0][dynamic_pressure]
    result_expr = result_pressure_expr.subs({density: density_, flow_speed: velocity_})
    return Quantity(result_expr)
