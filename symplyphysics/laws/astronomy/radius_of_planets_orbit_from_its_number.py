"""
Radius of planetary orbits from number
======================================

The Titius—Bode rule is an empirical formula that roughly describes the distances between the planets of the
Solar System and the Sun (the average radii of the orbits).

**Links:**

#. `Wikipedia, third formula <https://en.wikipedia.org/wiki/Titius%E2%80%93Bode_law#Original_formulation>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

orbit_radius = symbols.radius
"""
:symbols:`radius` of the planet's orbit.
"""

number_of_planet = symbols.whole_number
"""
Planet's number, starting from :math:`-1`. See :symbols:`whole_number`.
"""

first_constant = Quantity(0.4 * units.astronomical_unit, display_symbol="a")
"""
A quantity which is the free term in the formula.
"""

second_constant = Quantity(0.3 * units.astronomical_unit, display_symbol="b")
"""
A quantity which is The factor before the exponent.
"""

law = Eq(orbit_radius, first_constant + second_constant * 2**number_of_planet)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(number_of_planet_=number_of_planet)
@validate_output(orbit_radius)
def calculate_radius_of_orbit(number_of_planet_: int) -> Quantity:
    if number_of_planet_ < -1:
        raise ValueError("The planet number must be greater than -1 or equal.")

    result_expr = solve(law, orbit_radius, dict=True)[0][orbit_radius]
    result_expr = result_expr.subs({
        number_of_planet: number_of_planet_,
    })
    return Quantity(result_expr)
