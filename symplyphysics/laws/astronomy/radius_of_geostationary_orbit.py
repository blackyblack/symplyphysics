"""
Radius of geostationary orbit
=============================

A geostationary orbit is a circular orbit located above the Earth's equator (0° latitude),
where an artificial satellite orbits the planet with an angular velocity equal to the
angular speed of the Earth's rotation around its axis.

**Notation:**

#. :quantity_notation:`gravitational_constant`.

**Links:**

#. `Wikipedia, possible formula derivable from here <https://en.wikipedia.org/wiki/Geostationary_orbit#Derivation>`__.

..
    TODO: find link with exact formula
    TODO: move law to `gravity`?
"""

from sympy import (Eq, Rational, solve)
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
)

# Description

## Law is: R = (G * M / (w^2))^(1/3), where
## R - radius of geostationary orbit,
## G - gravitational constant,
## M - planet_mass of planet.
## w - angular speed of rotation of satellite.


orbital_radius = symbols.radius
"""
:symbols:`radius` of the satellite's geostationary orbit.
"""

planet_mass = symbols.mass
"""
:symbols:`mass`
"""

satellite_angular_speed = symbols.angular_speed
"""
:symbols:`angular_speed` of the satellite's rotation.
"""

law = Eq(orbital_radius,
    (quantities.gravitational_constant * planet_mass / (satellite_angular_speed**2))**Rational(1, 3))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_of_planet_=planet_mass, speed_rotation_satellite_=satellite_angular_speed)
@validate_output(orbital_radius)
def calculate_radius_of_orbit(mass_of_planet_: Quantity,
    speed_rotation_satellite_: Quantity) -> Quantity:
    result_expr = solve(law, orbital_radius, dict=True)[0][orbital_radius]
    result_expr = result_expr.subs({
        planet_mass: mass_of_planet_,
        satellite_angular_speed: speed_rotation_satellite_,
    })
    return Quantity(result_expr)
