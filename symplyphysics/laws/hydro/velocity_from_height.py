from sympy import (Eq, solve, sqrt)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input, validate_output)

# Description
## If hole appears in side wall or bottom of tank with liquid, liquid starts flowing out of this tank with some velocity.
## Velocity might be calculated with help of Torricelli's formula.
## Law: v = sqrt(2 * g * h), where
## v is velocity of liquid flowing out of the hole
## g is gravitation acceleration
## h is height of liquid above the hole.

# Conditions
## Liquid is ideal (no heat losses and no liquid friction losses)
## Hole is small (constant height at any point of the hole)

liquid_velocity = Symbol("liquid_velocity", units.velocity)
height_above_hole = Symbol("height_above_hole", units.length)

law = Eq(liquid_velocity, sqrt(2 * units.acceleration_due_to_gravity * height_above_hole))


def print() -> str:
    return print_expression(law)


@validate_input(height_=height_above_hole)
@validate_output(liquid_velocity)
def calculate_velocity(height_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, liquid_velocity, dict=True)[0][liquid_velocity]
    result_expr = result_velocity_expr.subs({height_above_hole: height_})
    return expr_to_quantity(result_expr)
