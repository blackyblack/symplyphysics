"""
Corrected planet period squared is proportional to cube of semimajor axis
=========================================================================

*Kepler's third law* of planetary motion relates the planet's *rotation period* to the *semi-major axis*
of its orbit. Kepler's laws use the assumption that the attracting mass is much greater than that of
the orbiting planet. If this assumption isn't made, the first and the second law stay the same, but
the third law requires a revision. As a result, the period of the planet's rotation depends not only
on the attracting mass, but also on the mass of the planet itself.
"""

from sympy import Eq, pi, solve
from sympy.physics.units import gravitational_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

rotation_period = Symbol("rotation_period", units.time)
"""
The planet's period of rotation.

Symbol:
    :code:`T`
"""

semimajor_axis = Symbol("semimajor_axis", units.length)
"""
The semi-major axis of the planet's orbit.

Symbol:
    :code:`a`
"""

attracting_mass = clone_symbol(symbols.basic.mass, display_symbol="M")
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the attracting body, such as the Sun.
"""

planetary_mass = clone_symbol(symbols.basic.mass)
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the orbiting planet.
"""

law = Eq(
    rotation_period**2,
    (4 * pi**2 * semimajor_axis**3) / (gravitational_constant * (attracting_mass + planetary_mass)),
)
r"""
:code:`T^2 = (4 * pi^2 * a^3) / (G * (M + m))`

Latex:
    .. math::
        T^2 = \frac{4 \pi^2 a^3}{G \left( M + m \right)}
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
