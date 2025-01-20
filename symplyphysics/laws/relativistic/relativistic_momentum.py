from sympy import (Eq, solve, sqrt)
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## Momentum (amount of motion) is a vector physical quantity that is a measure of the mechanical movement
## of a body. The relativistic momentum also takes into account speed limits equal to the speed of light.

## Law is: p = m * v / sqrt(1 - v^2 / c^2), where
## p - momentum,
## m - rest mass,
## v - velocity of the object in the selected frame of reference,
## c - speed of light in vacuum.

# Links: Wikipedia <https://en.wikipedia.org/wiki/Mass_in_special_relativity#Relativistic_energy%E2%80%93momentum_equation>

momentum = Symbol("momentum", units.momentum)

# This is equivalent of symbols.mass
rest_mass = Symbol("rest_mass", units.mass)
velocity = Symbol("velocity", units.velocity)

law = Eq(momentum, rest_mass * velocity / sqrt(1 - (velocity / speed_of_light)**2))


@validate_input(mass_=rest_mass, velocity_=velocity)
@validate_output(momentum)
def calculate_momentum(mass_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, momentum, dict=True)[0][momentum]
    result_expr = result_expr.subs({
        rest_mass: mass_,
        velocity: velocity_,
    })
    return Quantity(result_expr)
