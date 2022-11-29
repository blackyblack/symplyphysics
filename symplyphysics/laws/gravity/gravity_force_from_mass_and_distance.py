from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from sympy.physics.units import gravitational_constant

# Description
## Gravitational force between two object derivation from masses of objects and distance between them.
## Law: F = G * m1 * m2 / R**2, where
## F - gravitational force
## m1 and m2 - masses of objects
## R is distance between mass centers of objects
## G is gravitational constant

gravitational_force= symbols("gravitational_force")
mass1, mass2 = symbols("mass1, mass2")
distance_between_objects = symbols("distance")

law = Eq(gravitational_force, gravitational_constant * mass1 * mass2 / distance_between_objects**2)

def print():
    return pretty(law, use_unicode=False)

@validate_input(mass1_=units.mass, mass2_=units.mass, distance_between_objects_=units.length)
@validate_output(units.force)
def calculate_force_on_earth(mass1_: Quantity, mass2_: Quantity, distance_between_objects_: Quantity) -> Quantity:
    result_force_expr = solve(law, gravitational_force, dict=True)[0][gravitational_force]
    result_expr = result_force_expr.subs({mass1: mass1_, mass2: mass2_, distance_between_objects: distance_between_objects_})
    return expr_to_quantity(result_expr, 'gravity_force')
