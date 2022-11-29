from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from sympy.physics.units import gravitational_constant
from symplyphysics.laws.gravity import gravity_force_from_mass_and_distance as gravity_law
from symplyphysics.laws.dynamics import acceleration_from_force as newton2_law

# Description
## Every planet generates gravity field whichs causes free falling. Free fall acceleration depends on heigth above the planet surface.
## Law: g = G * M / (R+h)**2, where
## g is free fall acceleration
## G is gravitational constant
## M is mass of the planet
## R is radius of the planet
## h is heigth above the planet surface

# Conditions:
## Earth mass is 5.9742e24 kilos
## Earth radius is 6378.1e3 meters
## earth_mass = 5.9742e24 * units.kilogram
## earth_radius = 6378.1e3 * units.meter

free_fall_acceleration = symbols("free_fall_acceleration")
planet_mass = symbols("planet_mass")
planet_radius = symbols("planet_radius")
heigth_above_surface = symbols("heigth_above_surface")

## Let's define object's mass and discover how free fall acceleration depends on it
object_mass = symbols("object_mass")

law = Eq(free_fall_acceleration, gravitational_constant * planet_mass / (planet_radius + heigth_above_surface)**2)
gravitational_force = gravity_law.law.rhs.subs({
    gravity_law.mass1: planet_mass,
    gravity_law.mass2: object_mass,
    gravity_law.distance_between_objects: planet_radius + heigth_above_surface
    })

# This law might be easily derived from gravitational law via Newton's law #2
free_fall_acceleration = newton2_law.law.rhs.subs({
    newton2_law.force: gravitational_force,
    newton2_law.mass: object_mass
    })

def print():
    return pretty(law, use_unicode=False)
'''
@validate_input(generator_mass_=units.mass, generator_radius_=units.length)
@validate_output(units.acceleration)
def calculate_acceleration(generator_mass_: Quantity, generator_radius_: Quantity) -> Quantity:
    result_accel_expr = solve(law, acceleration_free_fall, dict=True)[0][acceleration_free_fall]
    result_expr = result_accel_expr.subs({generator_mass: generator_mass_, generator_radius: generator_radius_})
    return expr_to_quantity(result_expr, 'acceleration_free_fall')
    '''