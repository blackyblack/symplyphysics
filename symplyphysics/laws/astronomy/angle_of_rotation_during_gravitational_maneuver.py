r"""
Angle of rotation during gravitational maneuver
===============================================

A gravitational maneuver is a purposeful change in the trajectory and flight speed of a spacecraft under the influence
of the gravitational fields of celestial bodies.
The angle of the gravitational maneuver depends on the aiming range, the mass of the planet and the velocity of the rocket relative to the planet.

.. image:: https://upload.wikimedia.org/wikipedia/commons/a/ad/Gravity_assist\_-\_ru.svg
    :width: 240px

..
    TODO find link to law
    TODO move law to `gravity`?
"""

from sympy import Eq, solve, atan
from sympy.physics.units import gravitational_constant
from symplyphysics import (clone_as_symbol, symbols, Quantity,
    validate_input, validate_output)

angle = symbols.angle
"""
:symbols:`angle` of rotation during a gravitational maneuver (angle at which the velocity vector of the rocket rotates)
"""

planet_mass = clone_as_symbol(symbols.mass, display_symbol="M", display_latex="M")
"""
The :symbols:`mass` of the planet.
"""

aiming_range = clone_as_symbol(symbols.euclidean_distance, display_symbol="b", display_latex="b")
"""
The aiming range is the :symbols:`euclidean_distance` between the asymptote of the hyperbolic trajectory
of the circumnavigation of the planet and its focus coinciding with the center of the planet.
"""

rocket_speed = symbols.speed
"""
Rocket's :symbols:`speed` relative to the planet.
"""

law = Eq(angle, 2 * atan(gravitational_constant * planet_mass / (aiming_range * rocket_speed**2)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(planet_mass_=planet_mass, aiming_range_=aiming_range, rocket_speed_=rocket_speed)
@validate_output(angle)
def calculate_angle(planet_mass_: Quantity, aiming_range_: Quantity,
    rocket_speed_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, angle, dict=True)[0][angle]
    result_expr = result_velocity_expr.subs({
        planet_mass: planet_mass_,
        aiming_range: aiming_range_,
        rocket_speed: rocket_speed_
    })
    return Quantity(result_expr)
