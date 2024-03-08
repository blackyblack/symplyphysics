from sympy import Eq, solve, sqrt
from sympy.physics.units import gravitational_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.gravity import gravity_force_from_mass_and_distance as gravity_force_law
from symplyphysics.laws.dynamics import acceleration_from_force as acceleration_law
from symplyphysics.laws.kinematic import centripetal_acceleration_is_squared_velocity_by_radius as centripetal_law


# Law: V = √(G * M / (R + h))
# Where:
# V - initial velocity
# M - planet_mass of planet
# G - gravitational constant
# h - height above the planet surface
# R - radius of planet

velocity = Symbol("initial_velocity", units.velocity)
planet_mass = Symbol("planet_mass", units.mass)
radius = Symbol("radius", units.length)
height = Symbol("height", units.length)

law = Eq(velocity, sqrt(gravitational_constant * planet_mass / (radius + height)))

# This law might be derived via "gravity_force_from_mass_and_distance" law, "acceleration_from_force" law,
# "centripetal_acceleration_is_squared_velocity_by_radius" law.

# The radius of the orbit consists of the radius of the planet and the height above its surface.
centripetal_law_applied = centripetal_law.law.subs({
    centripetal_law.linear_velocity: velocity,
    centripetal_law.curve_radius: radius + height,
})
acceleration_derived = solve(centripetal_law_applied, centripetal_law.centripetal_acceleration, dict=True)[0][centripetal_law.centripetal_acceleration]

acceleration_law_applied = acceleration_law.law.subs({
    acceleration_law.acceleration: acceleration_derived,
})
force_derived = solve(acceleration_law_applied, acceleration_law.force, dict=True)[0][acceleration_law.force]

gravity_force_law_applied = gravity_force_law.law.subs({
    gravity_force_law.first_object_mass: planet_mass,
    gravity_force_law.gravitational_force: force_derived,
    gravity_force_law.second_object_mass: acceleration_law.mass,
    gravity_force_law.distance_between_mass_centers: radius + height,
})
velocity_derived = solve(gravity_force_law_applied, velocity, dict=True)[1][velocity]

# Check if derived velocity is same as declared.
assert expr_equals(velocity_derived, law.rhs)

def print_law() -> str:
    return print_expression(law)


@validate_input(planet_mass_=planet_mass, radius_=radius, height_=height)
@validate_output(velocity)
def calculate_velocity(planet_mass_: Quantity, radius_: Quantity, height_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, velocity, dict=True)[0][velocity]
    result_expr = result_velocity_expr.subs({
        planet_mass: planet_mass_,
        radius: radius_,
        height: height_
    })
    return Quantity(result_expr)
