from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from sympy.physics.units import gravitational_constant
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity

from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Every object generates gravity field around it. Any other object in this field is pulled toward generator.
## Gravitational force between two object if proportional to masses of objects and counter-proportional to distance between their mass centers.
## Law: F = G * m1 * m2 / R**2
## Where:
## F - gravitational force
## m1 and m2 - masses of objects
## R - distance between mass centers of objects
## G - gravitational constant

gravitational_force = Symbol("gravitational_force", units.force)
first_object_mass = Symbol("first_object_mass", units.mass)
second_object_mass = Symbol("second_object_mass", units.mass)
distance_between_mass_centers = Symbol("distance_between_mass_centers", units.length)

law = Eq(gravitational_force, gravitational_constant * first_object_mass * second_object_mass / distance_between_mass_centers**2)

def print(expr: Expr) -> str:
    symbols = [gravitational_force, first_object_mass, second_object_mass, distance_between_mass_centers]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(first_object_mass_=first_object_mass, second_object_mass_=second_object_mass, distance_between_objects_=distance_between_mass_centers)
@validate_output_symbol(gravitational_force)
def calculate_force(first_object_mass_: Quantity, second_object_mass_: Quantity, distance_between_objects_: Quantity) -> Quantity:
    result_force_expr = solve(law, gravitational_force, dict=True)[0][gravitational_force]
    result_expr = result_force_expr.subs({first_object_mass: first_object_mass_, second_object_mass: second_object_mass_, distance_between_mass_centers: distance_between_objects_})
    return expr_to_quantity(result_expr)
