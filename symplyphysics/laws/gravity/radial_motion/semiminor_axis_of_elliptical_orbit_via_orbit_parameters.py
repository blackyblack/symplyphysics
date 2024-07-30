"""
Semiminor axis of elliptical orbit via orbit parameters
=======================================================

The minor semiaxis can be found as a function of the sector speed of the planet,
the major semiaxis of its orbit, and the mass of body that attracts it, such as the Sun.
"""

from sympy import Eq, sqrt
from sympy.physics.units import gravitational_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    clone_symbol,
    symbols,
)

semiminor_axis = Symbol("semiminor_axis", units.length)
"""
The semi-minor axis of the planet's orbit.

Symbol:
    b
"""

sector_speed = Symbol("sector_speed", units.area / units.time)
r"""
The sector speed of the planet, i.e. the area swept by the planet per unit time.

Symbol:
    sigma

Latex:
    :math:`\sigma`
"""

semimajor_axis = Symbol("semimajor_axis", units.length)
"""
The semi-major axis of the planet's orbit.

Symbol:
    a
"""

attracting_mass = clone_symbol(symbols.basic.mass, "attracting_mass")
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the attracting body, such as the Sun.

Symbol:
    M
"""

law = Eq(
    semiminor_axis,
    2 * sector_speed * sqrt(semimajor_axis / (gravitational_constant * attracting_mass)),
)
r"""
b = 2 * sigma * sqrt(a / (G * M))

Latex:
    :math:`b = 2 \sigma \sqrt{\frac{a}{G M}}`
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
