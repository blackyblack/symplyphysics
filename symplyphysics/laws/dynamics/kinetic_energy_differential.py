from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The general formula for the kinetic energy of an object features its speed and momentum. This way it can be used
## not only in the case of variable mass, but also in the relativistic case.

# Law: dE = v * dp
## E - kinetic energy
## v - speed
## p - momentum (in general it is a function of speed)
## d(x) - exact differential of `x`

# Notes
## - Works in the case of special relativity as well

kinetic_energy_differential = Symbol("kinetic_energy_differential", units.energy)
speed = Symbol("speed", units.velocity)
momentum_differential = Symbol("momentum", units.momentum)

law = Eq(kinetic_energy_differential, speed * momentum_differential)

# TODO: derive from the differential definition of workand the generallized Newton's second law


@validate_input(
    speed_=speed,
    momentum_differential_=momentum_differential,
)
@validate_output(kinetic_energy_differential)
def calculate_kinetic_energy_differential(
    speed_: Quantity,
    momentum_differential_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        speed: speed_,
        momentum_differential: momentum_differential_,
    })
    return Quantity(result)
