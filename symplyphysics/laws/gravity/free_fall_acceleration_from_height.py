"""
Free fall acceleration from height
==================================

**Free fall** is any motion of a body where gravity is the only force acting upon it.
**Free fall acceleration** is the acceleration the body experiences during the free fall.

**Notation:**

#. :quantity_notation:`gravitational_constant`.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.gravity import gravity_force_from_mass_and_distance as gravity_law
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass as newton2_law

free_fall_acceleration = symbols.acceleration
"""
Free fall :symbols:`acceleration` of the body.
"""

planet_radius = symbols.radius
"""
:symbols:`radius` of the planet.
"""

elevation = symbols.height
"""
Elevation (:symbols:`height`) of the body from the planet's surface.
"""

planet_mass = symbols.mass
"""
:symbols:`mass` of the planet.
"""

law = Eq(free_fall_acceleration,
    quantities.gravitational_constant * planet_mass / (planet_radius + elevation)**2)
"""
:laws:symbol::

:laws:latex::
"""

# This law might be easily derived from gravitational law via Newton's law #2
## Distance between mass centers is radius of the planet plus height above it's surface.
_gravitational_force = gravity_law.law.rhs.subs({
    gravity_law.first_mass: planet_mass,
    gravity_law.distance_between_mass_centers: planet_radius + elevation
})

# Substitute mass first
_derived_free_fall_acceleration = newton2_law.law.rhs.subs(newton2_law.mass, gravity_law.second_mass)
_derived_free_fall_acceleration = _derived_free_fall_acceleration.subs(newton2_law.force,
    _gravitational_force)

# Check if derived acceleration is same as declared
assert expr_equals(_derived_free_fall_acceleration, law.rhs)


@validate_input(planet_mass_=planet_mass,
    planet_radius_=planet_radius,
    height_above_surface_=elevation)
@validate_output(free_fall_acceleration)
def calculate_acceleration(planet_mass_: Quantity, planet_radius_: Quantity,
    height_above_surface_: Quantity) -> Quantity:
    result_accel_expr = solve(law, free_fall_acceleration, dict=True)[0][free_fall_acceleration]
    result_expr = result_accel_expr.subs({
        planet_mass: planet_mass_,
        planet_radius: planet_radius_,
        elevation: height_above_surface_
    })
    return Quantity(result_expr)
