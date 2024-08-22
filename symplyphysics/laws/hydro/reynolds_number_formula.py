"""
Reynolds number formula
=======================

*Reynolds number* is a dimensionless quantity that characterizes the flow of a fluid.
It helps predict fluid flow patterns in different situations by measuring the ratio between inertial
and viscous forces. Low Reynolds numbers tend to correspond to laminar flows, while high Reynolds
numbers tend to correspond to turbulent flows.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

reynolds_number = Symbol("reynolds_number", dimensionless)
r"""
Reynolds number of the fluid.

Symbol:
    :code:`Re`

Latex:
    :math:`\text{Re}`
"""

characteristic_length = Symbol("characteristic_length", units.length)
"""
Characteristic length of the fluid container.

Symbol:
    :code:`L`
"""

density = Symbol("density", (units.mass / units.volume))
r"""
Density of the fluid.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

flow_speed = Symbol("flow_speed", units.velocity)
"""
Speed of the fluid flow.

Symbol:
    :code:`u`
"""

dynamic_viscosity = Symbol("dynamic_viscosity", units.pressure * units.time)
r"""
Dynamic viscosity of the fluid.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

law = Eq(reynolds_number, density * flow_speed * characteristic_length / dynamic_viscosity)
r"""
:code:`Re = rho * u * L / mu`

Latex:
    .. math::
        \text{Re} = \frac{\rho u L}{\mu}
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
    result = Quantity(result_applied)
    return convert_to_float(result)
