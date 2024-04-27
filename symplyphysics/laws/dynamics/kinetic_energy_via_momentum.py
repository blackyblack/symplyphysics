from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

# Description
## Kinetic energy can be expressed using linear momentum and mass.

# Law: E = p**2 / (2 * m)
## E - kinetic energy
## p - linear momentum
## m - mass

# Notes
## - This relation also holds in Quantum Mechanics.

kinetic_energy = Symbol("kinetic_energy", units.energy)
momentum = Symbol("momentum", units.momentum)
mass = symbols.basic.mass

law = Eq(kinetic_energy, momentum**2 / (2 * mass))

# TODO: derive law from kinetic energy and momentum expressions


@validate_input(
    momentum_=momentum,
    mass_=mass,
)
@validate_output(kinetic_energy)
def calculate_kinetic_energy(
    momentum_: Quantity,
    mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        momentum: momentum_,
        mass: mass_,
    })
    return Quantity(result)
