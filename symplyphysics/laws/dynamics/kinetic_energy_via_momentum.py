from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import (
    kinetic_energy_from_mass_and_velocity as energy_law,
)
from symplyphysics.definitions import (
    momentum_is_mass_times_velocity as momentum_def,
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

# Derive law from kinetic energy and momentum expressions

_energy_eqn = energy_law.law.subs({
    energy_law.kinetic_energy_of_body: kinetic_energy,
    energy_law.body_velocity: momentum_def.velocity,
    energy_law.mass: mass,
})

_momentum_eqn = momentum_def.definition.subs({
    momentum_def.momentum: momentum,
    momentum_def.mass: mass,
})

_kinetic_energy_expr = solve(
    (_energy_eqn, _momentum_eqn),
    (kinetic_energy, momentum_def.velocity),
    dict=True,
)[0][kinetic_energy]

assert expr_equals(_kinetic_energy_expr, law.rhs)


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