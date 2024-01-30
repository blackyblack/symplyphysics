from sympy import (Eq, solve, sqrt)
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Momentum (amount of motion) is a vector physical quantity that is a measure of the mechanical movement
## of a body. The relativistic angular momentum also takes into account speed limits equal to the speed of light.

## Law is: p = m * v / sqrt(1 - v^2 / c^2), where
## p - momentum,
## m - mass,
## v - velocity,
## c - speed of light in vacuum.

momentum = Symbol("momentum", units.momentum)

mass = Symbol("mass", units.mass)
velocity = Symbol("velocity", units.velocity)

law = Eq(momentum, mass * velocity / sqrt(1 - (velocity / speed_of_light)**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=mass, velocity_=velocity)
@validate_output(momentum)
def calculate_momentum(mass_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, momentum, dict=True)[0][momentum]
    result_expr = result_expr.subs({
        mass: mass_,
        velocity: velocity_,
    })
    return Quantity(result_expr)
