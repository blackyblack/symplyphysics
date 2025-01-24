"""
Speed of charged particles in gas via mobility
==============================================

In a gas, the free path length of charged particles is small and the particles, as they
move, experience many collisions in the volume between the electrodes. In this case, the
speed of directional movement will depend not on the potential difference traveled, but
on the local value of the electric intensity.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Electron_mobility#>`__.

..
    TODO: rename file
    TODO: create usual law `v = mu * E`
"""

from sympy import Eq, solve
from symplyphysics import (units, Quantity, SymbolNew, validate_input, validate_output, symbols)

speed = symbols.speed
"""
:symbols:`speed` of charged particles in the gas.
"""

mobility_at_unit_pressure = SymbolNew("mu", units.speed * units.pressure * units.length / units.voltage, display_latex="\\mu")
"""
Mobility of charged particles at unit pressure.
"""

pressure = symbols.pressure
"""
Pressure in the gas.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`
"""

law = Eq(speed, mobility_at_unit_pressure * electric_field_strength / pressure)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mobility_at_unit_pressure_=mobility_at_unit_pressure,
    pressure_=pressure,
    electric_intensity_=electric_field_strength)
@validate_output(speed)
def calculate_velocity(mobility_at_unit_pressure_: Quantity, pressure_: Quantity,
    electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, speed, dict=True)[0][speed]
    result_expr = result_expr.subs({
        mobility_at_unit_pressure: mobility_at_unit_pressure_,
        pressure: pressure_,
        electric_field_strength: electric_intensity_,
    })
    return Quantity(result_expr)
