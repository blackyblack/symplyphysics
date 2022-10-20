from turtle import distance
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from sympy.physics.units import gravitational_constant

# Description
## Field strength of gravity field E = G * M / R^2, where
## G is the gravitational constant
## M is the mass of field generator
## R is the distance to the generator's mass center
## Gravity force applied to object F = E * m, where
## E is the gravity field strength
## m is mass of the object.
## So gravity field strength in any point is gravity acceleration in this point.

gravity_force, generator_mass, object_mass, distance_between_mass_centers  = symbols('gravity_force generator_mass object_mass distance')
law = Eq(gravity_force, gravitational_constant * generator_mass * object_mass / (distance_between_mass_centers ** 2))

def print():
    return pretty(law, use_unicode=False)

@validate_input(generator_mass_=units.mass, object_mass_=units.mass, distance_between_mass_centers_=units.length)
@validate_output(units.force)
def calculate_force(generator_mass_: Quantity, object_mass_: Quantity, distance_between_mass_centers_: Quantity) -> Quantity:
    result_force_expr = solve(law, gravity_force, dict=True)[0][gravity_force]
    result_expr = result_force_expr.subs({generator_mass: generator_mass_, object_mass: object_mass_, distance_between_mass_centers: distance_between_mass_centers_})
    return expr_to_quantity(result_expr, 'gravity_force')
