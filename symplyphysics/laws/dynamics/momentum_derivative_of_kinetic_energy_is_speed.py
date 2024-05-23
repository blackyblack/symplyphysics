from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)

# Description
## The general formula for the kinetic energy of an object features its speed and momentum. This way it can be used
## not only in the case of variable mass, but also in the relativistic case.

# Law: dE(p(v))/dp(v) = v
## E - kinetic energy
## p - momentum
## v - speed
## d/dp - derivative w.r.t. momentum

# Notes
## - Works in the case of special relativity as well

kinetic_energy = Function("kinetic_energy", units.energy)
momentum = Function("momentum", units.momentum)
speed = Symbol("speed", units.velocity)

law = Eq(
    Derivative(kinetic_energy(momentum(speed)), momentum(speed)),
    speed,
)

# TODO: derive from the differential definition of work and the generalized Newton's second law


@validate_input(
    kinetic_energy_change_=kinetic_energy,
    momentum_change_=momentum,
)
@validate_output(speed)
def calculate_speed(
    kinetic_energy_change_: Quantity,
    momentum_change_: Quantity,
) -> Quantity:
    kinetic_energy_ = (kinetic_energy_change_ / momentum_change_) * momentum(speed)
    result = law.lhs.subs(
        kinetic_energy(momentum(speed)),
        kinetic_energy_,
    ).doit()
    return Quantity(result)
