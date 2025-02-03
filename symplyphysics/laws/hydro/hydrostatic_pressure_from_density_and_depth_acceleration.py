"""
Hydrostatic pressure from density, height and acceleration
==========================================================

Hydrostatic pressure is the pressure exerted by a fluid at equilibrium at a given point
within the fluid, due to the force of gravity.

**Conditions:**

#. The only force acting on the fluid is the gravitational force.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import symbols, Quantity, validate_input, validate_output

density = symbols.density
"""
:symbols:`density` of the fluid.
"""

height = symbols.height
"""
:symbols:`height` of the fluid column.
"""

acceleration = symbols.acceleration
"""
:symbols:`acceleration` of the vessel.
"""

hydrostatic_pressure = symbols.hydrostatic_pressure
"""
:symbols:`hydrostatic_pressure` of the fluid.
"""

law = Eq(hydrostatic_pressure, density * acceleration * height)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(density_=density, depth_=height, acceleration_=acceleration)
@validate_output(hydrostatic_pressure)
def calculate_hydrostatic_pressure(density_: Quantity, depth_: Quantity,
    acceleration_: Quantity) -> Quantity:
    result_pressure_expr = solve(law, hydrostatic_pressure, dict=True)[0][hydrostatic_pressure]
    result_expr = result_pressure_expr.subs({
        density: density_,
        height: depth_,
        acceleration: acceleration_,
    })
    return Quantity(result_expr)
