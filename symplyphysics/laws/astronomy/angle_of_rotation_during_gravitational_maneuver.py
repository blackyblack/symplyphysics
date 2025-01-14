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
from symplyphysics import symbols, Quantity, validate_input, validate_output, quantities

angle = symbols.angle
"""
:symbols:`angle` of rotation during a gravitational maneuver (angle at which the velocity vector of the rocket rotates)
"""

planet_mass = symbols.mass
"""
The :symbols:`mass` of the planet.
"""

aiming_range = symbols.euclidean_distance
"""
The aiming range is the :symbols:`euclidean_distance` between the asymptote of the hyperbolic trajectory
of the circumnavigation of the planet and its focus coinciding with the center of the planet.
"""

rocket_speed = symbols.speed
"""
Rocket's :symbols:`speed` relative to the planet.
"""

law = Eq(angle, 2 * atan(quantities.gravitational_constant * planet_mass / (aiming_range * rocket_speed**2)))
"""
:code:`phi = 2 * atan(G * m / (d * v^2))`

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
