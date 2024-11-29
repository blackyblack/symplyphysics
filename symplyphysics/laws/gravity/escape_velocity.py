"""
First escape speed
==================

**First escape speed** is the minimum speed needed for an object to escape from contact with or orbit of a primary body.

**Notation:**

#. :quantity_notation:`gravitational_constant`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Escape_velocity#From_an_orbiting_body>`__.

..
    TODO rename file
"""

from sympy import Eq, solve, sqrt
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.gravity import gravity_force_from_mass_and_distance as gravity_force_law
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass as acceleration_law
from symplyphysics.laws.kinematics import centripetal_acceleration_via_linear_speed_and_radius as centripetal_law

speed = symbols.speed
"""
Escape :symbols:`speed` of the body.
"""

radius = symbols.radius
"""
:symbols:`radius` of the planet.
"""

height = symbols.height
"""
Elevation (:symbols:`height`) of the body from the surface of the planet
"""

planet_mass = symbols.mass
"""
:symbols:`mass` of the planet.
"""

law = Eq(speed, sqrt(quantities.gravitational_constant * planet_mass / (radius + height)))
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via "gravity_force_from_mass_and_distance" law, "acceleration_from_force" law,
# "centripetal_acceleration_is_squared_velocity_by_radius" law.

# The radius of the orbit consists of the radius of the planet and the height above its surface.
_centripetal_law_applied = centripetal_law.law.subs({
    centripetal_law.speed: speed,
    centripetal_law.radius_of_curvature: radius + height,
})
_acceleration_derived = solve(_centripetal_law_applied,
    centripetal_law.centripetal_acceleration,
    dict=True)[0][centripetal_law.centripetal_acceleration]

_acceleration_law_applied = acceleration_law.law.subs({
    acceleration_law.symbols.acceleration: _acceleration_derived,
})
_force_derived = solve(_acceleration_law_applied, acceleration_law.symbols.force,
    dict=True)[0][acceleration_law.symbols.force]

# Let's write down Newton's second law: ma = F. F is, in this case, the force of gravity. And in the general case,
# when a body moves along a circle with a constant speed in modulus, its acceleration is equal to the centripetal
# acceleration.
_gravity_force_law_applied = gravity_force_law.law.subs({
    gravity_force_law.first_mass: planet_mass,
    gravity_force_law.gravitational_force: _force_derived,
    gravity_force_law.second_mass: acceleration_law.symbols.mass,
    gravity_force_law.distance_between_mass_centers: radius + height,
})
# The first cosmic speed is the minimum horizontal speed that must be given to an object so that it moves in
# a circular orbit around the planet. Based on this definition, the speed could be left negative, implying a different
# direction. But most often they are written with a plus sign, implying the modulus of the horizontal component of the
# speed tangent to the orbit. Therefore, the first solution that returns a minus is ignored.
_velocity_derived = solve(_gravity_force_law_applied, speed, dict=True)[1][speed]

# Check if derived speed is same as declared.
assert expr_equals(_velocity_derived, law.rhs)


@validate_input(planet_mass_=planet_mass, radius_=radius, height_=height)
@validate_output(speed)
def calculate_velocity(planet_mass_: Quantity, radius_: Quantity, height_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, speed, dict=True)[0][speed]
    result_expr = result_velocity_expr.subs({
        planet_mass: planet_mass_,
        radius: radius_,
        height: height_
    })
    return Quantity(result_expr)
