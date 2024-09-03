"""
Dynamic pressure from speed
===========================

*Dynamic pressure*, sometimes called *velocity pressure*, is a physical quantity denoting the
pressure caused by a flowing fluid. It is numerically equal to the kinetic energy of the
fluid per unit volume (see :doc:`laws.quantities.quantity_is_volumetric_density_times_volume`).

**Notes:**

#. Many authors define this quantity only for *incompressible flows*, but others extend it for
   compressible flows as well.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)

dynamic_pressure = Symbol("dynamic_pressure", units.pressure)
"""
Dynamic pressure of the fluid.

Symbol:
    :code:`q`
"""

density = Symbol("density", units.mass / units.volume)
r"""
Density of the fluid.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

flow_speed = Symbol("flow_speed", units.velocity)
"""
Flow speed of the fluid.

Symbol:
    :code:`u`
"""

law = Eq(dynamic_pressure, density * flow_speed**2 / 2)
r"""
:code:`q = 1/2 * rho * u^2`

Latex:
    .. math::
        q = \frac{1}{2} \rho u^2
"""


@validate_input(density_=density, velocity_=flow_speed)
@validate_output(dynamic_pressure)
def calculate_pressure(density_: Quantity, velocity_: Quantity) -> Quantity:
    result_pressure_expr = solve(law, dynamic_pressure, dict=True)[0][dynamic_pressure]
    result_expr = result_pressure_expr.subs({density: density_, flow_speed: velocity_})
    return Quantity(result_expr)
