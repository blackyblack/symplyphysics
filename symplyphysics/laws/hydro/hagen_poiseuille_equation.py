"""
Hagen—Poiseuille equation
========================

In non-ideal fluid mechanics, the Hagen—Poiseuille equation is a physical law that gives
the pressure drop in an incompressible and Newtonian fluid in laminar flow flowing
through a long cylindrical pipe of constant cross section. It can be successfully
applied to air flow in lung alveoli, or the flow through a drinking straw.

**Conditions:**

#. The fluid is incompressible and Newtonian.
#. The flow of the fluid is laminar.
#. The pipe is long enough for the flow to be laminar.
#. The pipe has constant cross section.

**Links:**

#. `Wikipedia, first part of the first equation <https://en.wikipedia.org/wiki/Hagen%E2%80%93Poiseuille_equation#Equation>`__.
"""

from sympy import Eq, solve, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

dynamic_viscosity = symbols.dynamic_viscosity
"""
:symbols:`dynamic_viscosity` of the fluid.
"""

length = symbols.length
"""
:symbols:`length` of the pipe.
"""

flow_rate = symbols.volumetric_flow_rate
"""
:symbols:`volumetric_flow_rate` of the fluid through the pipe.
"""

radius = symbols.radius
"""
:symbols:`radius` of the pipe.
"""

pressure_change = clone_as_symbol(symbols.pressure, display_symbol="Delta(p)", display_latex="\\Delta p")
"""
Change in :symbols:`pressure` between the two ends.
"""

law = Eq(pressure_change, 8 * dynamic_viscosity * length * flow_rate / (pi * radius**4))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    dynamic_viscosity_=dynamic_viscosity,
    length_=length,
    flow_rate_=flow_rate,
    radius_=radius,
)
@validate_output(pressure_change)
def calculate_delta_pressure(dynamic_viscosity_: Quantity, length_: Quantity, flow_rate_: Quantity,
    radius_: Quantity) -> Quantity:
    result_expr = solve(law, pressure_change, dict=True)[0][pressure_change]
    result_applied = result_expr.subs({
        dynamic_viscosity: dynamic_viscosity_,
        length: length_,
        flow_rate: flow_rate_,
        radius: radius_,
    })
    return Quantity(result_applied)
