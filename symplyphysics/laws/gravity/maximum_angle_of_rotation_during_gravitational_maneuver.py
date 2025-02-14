"""
Maximum angle of rotation during gravitational maneuver
=======================================================

A gravitational maneuver is a purposeful change in the trajectory and flight speed of a
spacecraft under the influence of the gravitational fields of celestial bodies. The
maximum angle of rotation of the rocket during a gravitational maneuver depends on the
first cosmic velocity of the planet around which the maneuver is performed, and the speed
of the rocket relative to this planet. The law describes the closest possible maneuver,
when aiming range is at planet radius.

**Conditions:**

#. Aiming range equals the planet radius.

..
    TODO find link
"""

from sympy import (Eq, solve, atan)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

maximum_angle = symbols.angle
"""
Maximum :symbols:`angle` of rotation during the gravitational maneuver.
"""

first_cosmic_speed = clone_as_symbol(symbols.speed, subscript="1")
"""
First cosmic :symbols:`speed` of the planet.
"""

rocket_speed = symbols.speed
"""
:symbols:`speed` relative to the planet.
"""

law = Eq(maximum_angle, atan((first_cosmic_speed / rocket_speed)**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(first_cosmic_velocity_planet_=first_cosmic_speed, rocket_speed_=rocket_speed)
@validate_output(maximum_angle)
def calculate_maximum_angle(first_cosmic_velocity_planet_: Quantity,
    rocket_speed_: Quantity) -> Quantity:
    result_expr = solve(law, maximum_angle, dict=True)[0][maximum_angle]
    result_expr = result_expr.subs({
        first_cosmic_speed: first_cosmic_velocity_planet_,
        rocket_speed: rocket_speed_,
    })
    return Quantity(result_expr)
