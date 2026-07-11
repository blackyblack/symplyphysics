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
"""

from sympy import (Eq, Rational, solve, Symbol as SymSymbol)
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.classical_mechanics.dynamics.force import acceleration_is_force_over_mass as newtons_second_law
from symplyphysics.classical_mechanics.dynamics.gravity import gravity_force_from_mass_and_distance as gravity_law
from symplyphysics.classical_mechanics.kinematics.rotational_motion import centripetal_acceleration_via_angular_speed_and_radius as centripetal_law

orbital_radius = symbols.radius
"""
:symbols:`radius` of the satellite's geostationary orbit.
"""

planet_mass = symbols.mass
"""
:symbols:`mass` of the attracting body (planet).
"""

satellite_angular_speed = symbols.angular_speed
"""
:symbols:`angular_speed` of the satellite's rotation.
"""

law = Eq(orbital_radius,
    (quantities.gravitational_constant * planet_mass / (satellite_angular_speed**2))**Rational(
    1, 3))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law from Newton's second law of motion applied to the circular motion of the
# satellite in the gravitational field of the planet.

_radius = SymSymbol("radius", positive=True)
_satellite_mass = newtons_second_law.mass

## The satellite moves in a circular orbit, so its acceleration is centripetal.
_acceleration_expr = centripetal_law.law.rhs.subs({
    centripetal_law.angular_speed: satellite_angular_speed,
    centripetal_law.radius_of_curvature: _radius,
})

## The only force acting on the satellite is the gravitational pull of the planet.
_force_expr = gravity_law.law.rhs.subs({
    gravity_law.first_mass: planet_mass,
    gravity_law.second_mass: _satellite_mass,
    gravity_law.distance_between_mass_centers: _radius,
})

_newtons_eqn = newtons_second_law.law.subs({
    newtons_second_law.acceleration: _acceleration_expr,
    newtons_second_law.force: _force_expr,
    newtons_second_law.mass: _satellite_mass,
})

_radius_derived = solve(_newtons_eqn, _radius)[0]

assert expr_equals(_radius_derived, law.rhs)


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
