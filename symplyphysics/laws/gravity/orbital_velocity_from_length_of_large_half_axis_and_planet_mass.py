"""
Orbital speed from semimajor axis and planet mass
=================================================

The orbital speed of a body is the speed at which it rotates around the barycenter of the
system, usually around a more massive body. It can be calculated from the mass of the
planet and the orbit configuration parameters.

**Notation:**

#. :quantity_notation:`gravitational_constant`.

**Conditions:**

#. The distance to the barycenter :math:`r` must not exceed the length :math:`a` of the
   semi-major axis of the orbit.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Orbital_speed#Instantaneous_orbital_speed>`__.
"""

from sympy import Eq, solve, sqrt
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities

orbital_speed = symbols.speed
"""
Orbital :symbols:`speed` of the satellite.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` at which the speed is calculated.
"""

semimajor_axis = symbols.semimajor_axis
"""
:symbols:`semimajor_axis` of the satellite's orbit.
"""

planet_mass = symbols.mass
"""
:symbols:`mass` of the planet.
"""

law = Eq(
    orbital_speed,
    sqrt(quantities.gravitational_constant * planet_mass * ((2 / distance) - (1 / semimajor_axis))))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(planet_mass_=planet_mass,
    distance_=distance,
    large_half_axis_length_=semimajor_axis)
@validate_output(orbital_speed)
def calculate_orbital_velocity(planet_mass_: Quantity, distance_: Quantity,
    large_half_axis_length_: Quantity) -> Quantity:
    if distance_.scale_factor > large_half_axis_length_.scale_factor:
        raise ValueError(
            "The distance between the rotating body and the central body must be less or equal to the large half-axis."
        )
    result_velocity_expr = solve(law, orbital_speed, dict=True)[0][orbital_speed]
    result_expr = result_velocity_expr.subs({
        planet_mass: planet_mass_,
        distance: distance_,
        semimajor_axis: large_half_axis_length_
    })
    return Quantity(result_expr)
