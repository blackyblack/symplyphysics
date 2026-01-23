"""
Semiminor axis of elliptical orbit via orbit parameters
=======================================================

The minor semiaxis can be found as a function of the sector speed of the planet,
the major semiaxis of its orbit, and the mass of body that attracts it, such as the Sun.

**Notation:**

#. :quantity_notation:`gravitational_constant`.

**Links:**

#. Sivukhin, D.V. (1979). *Obshchiy kurs fiziki* [General course of Physics], vol. 1, p. 318, (58.5)
"""

from sympy import Eq, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
    symbols,
    quantities,
)

semiminor_axis = symbols.semiminor_axis
"""
The :symbols:`semiminor_axis` of the planet's orbit.
"""

sector_speed = symbols.sector_speed
"""
The :symbols:`sector_speed` of the planet, i.e. the area swept by the planet per unit time.
"""

semimajor_axis = symbols.semimajor_axis
"""
The :symbols:`semimajor_axis` of the planet's orbit.
"""

attracting_mass = clone_as_symbol(symbols.mass, display_symbol="M", display_latex="M")
"""
The :symbols:`mass` of the attracting body, such as the Sun.
"""

law = Eq(
    semiminor_axis,
    2 * sector_speed * sqrt(semimajor_axis / (quantities.gravitational_constant * attracting_mass)),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    sector_speed_=sector_speed,
    semimajor_axis_=semimajor_axis,
    attracting_body_mass_=attracting_mass,
)
@validate_output(semiminor_axis)
def calculate_semiminor_axis(
    sector_speed_: Quantity,
    semimajor_axis_: Quantity,
    attracting_body_mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        sector_speed: sector_speed_,
        semimajor_axis: semimajor_axis_,
        attracting_mass: attracting_body_mass_,
    })
    return Quantity(result)
