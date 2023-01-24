from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity, SI, sqrt
)
from sympy.physics.units import acceleration_due_to_gravity as g

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

liquid_velocity, height_above_hole = symbols('liquid_velocity height_above_hole')
law = Eq(liquid_velocity, sqrt(2 * units.acceleration_due_to_gravity * height_above_hole))

def print():
    return pretty(law, use_unicode=False)

@validate_input(height_ = units.length)
@validate_output(units.velocity)
def calculate_velocity(height_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, liquid_velocity, dict=True)[0][liquid_velocity]        
    result_expr = result_velocity_expr.subs({height_above_hole: height_})
    return expr_to_quantity(result_expr, 'liquid_velocity')
