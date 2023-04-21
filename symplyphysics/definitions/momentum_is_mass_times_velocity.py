from sympy import (Eq, solve)
from symplyphysics import (
    units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol
)

# Description
## Momentum is the multiplication of velocity and mass. As velocity is vector, momentum is vector as well and it is collinear with velocity.

# Definition: P = m * V
# Where:
## m is a mass of the object
## V is it's velocity
## P is momentum

momentum = Symbol("momentum", units.momentum)
mass = Symbol("mass", units.mass)
velocity = Symbol("velocity", units.velocity)

definition = Eq(momentum, mass * velocity)

definition_units_SI = units.kilogram * units.meter / units.second

def print() -> str:
    return print_expression(definition)

@validate_input_symbols(velocity_=velocity, mass_=mass)
@validate_output_symbol(momentum)
def calculate_momentum(mass_: Quantity, velocity_: Quantity) -> Quantity:
    solved = solve(definition, momentum, dict=True)[0][momentum]
    result_expr = solved.subs({mass: mass_, velocity: velocity_})
    return expr_to_quantity(result_expr)
