from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from sympy.physics.units import gravitational_constant

# Description
## Every object generates gravity field around it. Any other object in this field is pulled toward generator.
## Gravitational force between two object if proportional to masses of objects and counter-proportional to distance between their mass centers.
## Law: F = G * m1 * m2 / R**2
## Where:
## F - gravitational force
## m1 and m2 - masses of objects
## R - distance between mass centers of objects
## G - gravitational constant

gravitational_force= symbols("gravitational_force")
first_object_mass, second_object_mass = symbols("first_object_mass second_object_mass")
distance_between_mass_centers = symbols("distance")

law = Eq(gravitational_force, gravitational_constant * first_object_mass * second_object_mass / distance_between_mass_centers**2)

def print():
    return pretty(law, use_unicode=False)

@validate_input(first_object_mass_=units.mass, second_object_mass_=units.mass, distance_between_objects_=units.length)
@validate_output(units.force)
def calculate_force(first_object_mass_: Quantity, second_object_mass_: Quantity, distance_between_objects_: Quantity) -> Quantity:
    result_force_expr = solve(law, gravitational_force, dict=True)[0][gravitational_force]
    result_expr = result_force_expr.subs({first_object_mass: first_object_mass_, second_object_mass: second_object_mass_, distance_between_mass_centers: distance_between_objects_})
    return expr_to_quantity(result_expr, "gravity_force")
