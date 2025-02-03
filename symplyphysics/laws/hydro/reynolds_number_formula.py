"""
Reynolds number formula
=======================

*Reynolds number* is a dimensionless quantity that characterizes the flow of a fluid. It
helps predict fluid flow patterns in different situations by measuring the ratio between
inertial and viscous forces. Low Reynolds numbers tend to correspond to laminar flows,
while high Reynolds numbers tend to correspond to turbulent flows.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Reynolds_number#Definition>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
)

reynolds_number = symbols.reynolds_number
"""
:symbols:`reynolds_number` of the fluid.
"""

characteristic_length = symbols.characteristic_length
"""
:symbols:`characteristic_length` of the system.
"""

density = symbols.density
"""
:symbols:`density` of the fluid.
"""

flow_speed = symbols.flow_speed
"""
:symbols:`flow_speed`
"""

dynamic_viscosity = symbols.dynamic_viscosity
"""
:symbols:`dynamic_viscosity` of the fluid.
"""

law = Eq(reynolds_number, density * flow_speed * characteristic_length / dynamic_viscosity)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(diameter_=characteristic_length,
    density_=density,
    velocity_=flow_speed,
    dynamic_viscosity_=dynamic_viscosity)
@validate_output(reynolds_number)
def calculate_reynolds_number(diameter_: Quantity, density_: Quantity, velocity_: Quantity,
    dynamic_viscosity_: Quantity) -> float:
    result_expr = solve(law, reynolds_number, dict=True)[0][reynolds_number]
    result_applied = result_expr.subs({
        characteristic_length: diameter_,
        density: density_,
        flow_speed: velocity_,
        dynamic_viscosity: dynamic_viscosity_
    })
    return convert_to_float(result_applied)
