"""
Corrected planet period squared is proportional to cube of semimajor axis
=========================================================================

*Kepler's third law* of planetary motion relates the planet's *rotation period* to the *semi-major axis*
of its orbit. Kepler's laws use the assumption that the attracting mass is much greater than that of
the orbiting planet. If this assumption isn't made, the first and the second law stay the same, but
the third law requires a revision. As a result, the period of the planet's rotation depends not only
on the attracting mass, but also on the mass of the planet itself.

**Notation:**

#. :quantity_notation:`gravitational_constant`.
"""

from sympy import Eq, pi, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)

rotation_period = symbols.period
"""
The planet's :symbols:`period` of rotation.
"""

semimajor_axis = symbols.semimajor_axis
"""
The :symbols:`semimajor_axis` of the planet's orbit.
"""

attracting_mass = clone_as_symbol(symbols.mass, display_symbol="M", display_latex="M")
"""
The :symbols:`mass` of the attracting body, such as the Sun.
"""

planetary_mass = symbols.mass
"""
The :symbols:`mass` of the orbiting planet.
"""

law = Eq(
    rotation_period**2,
    (4 * pi**2 * semimajor_axis**3) / (quantities.gravitational_constant * (attracting_mass + planetary_mass)),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    semimajor_axis_=semimajor_axis,
    attracting_mass_=attracting_mass,
    planetary_mass_=planetary_mass,
)
@validate_output(rotation_period)
def calculate_rotation_period(
    semimajor_axis_: Quantity,
    attracting_mass_: Quantity,
    planetary_mass_: Quantity,
) -> Quantity:
    expr = solve(law, rotation_period)[1]
    result = expr.subs({
        semimajor_axis: semimajor_axis_,
        attracting_mass: attracting_mass_,
        planetary_mass: planetary_mass_,
    })
    return Quantity(result)
