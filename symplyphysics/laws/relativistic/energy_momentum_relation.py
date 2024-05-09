from sympy import Eq, solve
from symplyphysics import (
    symbols,
    clone_symbol,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The energy—momentum relation, also called relativistic dispersion relation, is a relativistic
## equation relating total energy to invariant mass and momentum. It is the extension of
## [mass-energy equivalence](./energy_is_mass.py) for bodies or systems with non-zero momentum.

# Law: E**2 = (p * c)**2 + (m0 * c**2)**2
## E - total, or relativistic, energy
## p - momentum
## m0 - rest, or invariant, mass
## c - speed of light

relativistic_energy = Symbol("relativistic_energy", units.energy)
relativistic_momentum = Symbol("relativistic_momentum", units.momentum)
invariant_mass = clone_symbol(symbols.basic.mass, "invariant_mass")

law = Eq(
    relativistic_energy**2,
    (relativistic_momentum * units.speed_of_light)**2
    + (invariant_mass * units.speed_of_light**2)**2,
)


@validate_input(
    relativistic_momentum_=relativistic_momentum,
    invariant_mass_=invariant_mass,
)
@validate_output(relativistic_energy)
def calculate_relativistic_energy(
    relativistic_momentum_: Quantity,
    invariant_mass_: Quantity,
) -> Quantity:
    expr = solve(law, relativistic_energy)[1]
    result = expr.subs({
        relativistic_momentum: relativistic_momentum_,
        invariant_mass: invariant_mass_,
    })
    return Quantity(result)
