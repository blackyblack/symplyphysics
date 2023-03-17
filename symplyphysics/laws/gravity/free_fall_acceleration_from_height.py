from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units, simplify,
    validate_input, validate_output, expr_to_quantity
)
from sympy.physics.units import gravitational_constant
from symplyphysics.laws.gravity import gravity_force_from_mass_and_distance as gravity_law
from symplyphysics.laws.dynamics import acceleration_from_force as newton2_law

# Description
## Every planet generates gravity field which causes free falling. Free fall acceleration depends on height above the planet surface.
## Law: g = G * M / (R + h)**2
## Where:
## g is free fall acceleration
## G is gravitational constant
## M is mass of the planet
## R is radius of the planet
## h is height above the planet surface

free_fall_acceleration = symbols("free_fall_acceleration")
planet_mass = symbols("planet_mass")
planet_radius = symbols("planet_radius")
height_above_surface = symbols("height_above_surface")

law = Eq(free_fall_acceleration, gravitational_constant * planet_mass / (planet_radius + height_above_surface)**2)

# This law might be easily derived from gravitational law via Newton's law #2
## Distance between mass centers is radius of the planet plus height above it's surface.
gravitational_force = gravity_law.law.rhs.subs({
    gravity_law.first_object_mass: planet_mass,    
    gravity_law.distance_between_mass_centers: planet_radius + height_above_surface
    })

derived_free_fall_acceleration = newton2_law.law.rhs.subs({
    newton2_law.force: gravitational_force,
    newton2_law.mass: gravity_law.second_object_mass
    })

# Check if derived acceleration is same as declared    
difference = simplify(derived_free_fall_acceleration - law.rhs)
assert(difference == 0)

def print():
    return pretty(law, use_unicode=False)
    
@validate_input(planet_mass_=units.mass, planet_radius_=units.length, height_above_surface_=units.length)
@validate_output(units.acceleration)
def calculate_acceleration(planet_mass_: Quantity, planet_radius_: Quantity, height_above_surface_: Quantity) -> Quantity:
    result_accel_expr = solve(law, free_fall_acceleration, dict=True)[0][free_fall_acceleration]
    result_expr = result_accel_expr.subs({planet_mass: planet_mass_, planet_radius: planet_radius_, height_above_surface: height_above_surface_})
    return expr_to_quantity(result_expr, "free_fall_acceleration")
