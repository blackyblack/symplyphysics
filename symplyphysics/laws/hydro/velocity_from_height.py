from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity, sqrt
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

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

def print(expr: Expr) -> str:
    symbols = [liquid_velocity, height_above_hole]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(height_=height_above_hole)
@validate_output_symbol(liquid_velocity)
def calculate_velocity(height_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, liquid_velocity, dict=True)[0][liquid_velocity]        
    result_expr = result_velocity_expr.subs({height_above_hole: height_})
    return expr_to_quantity(result_expr)
