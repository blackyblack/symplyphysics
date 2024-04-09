from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols)

# Description
## Momentum is the multiplication of velocity and mass. As velocity is vector, momentum is vector as well and it is collinear with velocity.

# Definition: P = m * V
# Where:
## m is a mass of the object
## V is it's velocity
## P is momentum

momentum = Symbol("momentum", units.momentum)
velocity = Symbol("velocity", units.velocity)

definition = Eq(momentum, symbols.basic.mass * velocity)

definition_units_SI = units.kilogram * units.meter / units.second


def print_law() -> str:
    return print_expression(definition)


@validate_input(velocity_=velocity, mass_=symbols.basic.mass)
@validate_output(momentum)
def calculate_momentum(mass_: Quantity, velocity_: Quantity) -> Quantity:
    solved = solve(definition, momentum, dict=True)[0][momentum]
    result_expr = solved.subs({symbols.basic.mass: mass_, velocity: velocity_})
    return Quantity(result_expr)
