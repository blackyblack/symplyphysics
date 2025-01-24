"""
Radius of planetary orbits from number
======================================

The Titius—Bode rule is an empirical formula that roughly describes the distances between the planets of the
Solar System and the Sun (the average radii of the orbits).

**Notes:**

#. This rule fails to predict the correct value of the Neptune's orbit radius.

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

planet_number = symbols.whole_number
"""
The value :math:`-\\infty` corresponts to Mercury, :math:`0` to Venus, :math:`1` to
Earth, :math:`2` to Mars, :math:`3` to Ceres, :math:`4` to Jupyter, :math:`5` to Saturn,
:math:`6` to Uranus, and :math:`7` to Pluto.
"""

first_constant = Quantity(0.4 * units.astronomical_unit, display_symbol="a")
"""
A quantity which is the free term in the formula.
"""

second_constant = Quantity(0.3 * units.astronomical_unit, display_symbol="b")
"""
A quantity which is the factor before the exponent.
"""

law = Eq(orbit_radius, first_constant + second_constant * 2**planet_number)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(planet_number_=planet_number)
@validate_output(orbit_radius)
def calculate_radius_of_orbit(planet_number_: int) -> Quantity:
    if planet_number_ < -1:
        raise ValueError("The planet number must be greater than -1 or equal.")

    result_expr = solve(law, orbit_radius, dict=True)[0][orbit_radius]
    result_expr = result_expr.subs({
        planet_number: planet_number_,
    })
    return Quantity(result_expr)
