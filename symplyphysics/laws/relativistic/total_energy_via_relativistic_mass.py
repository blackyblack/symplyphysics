from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

# Description
## Fundamentally the energy of an object is synonimical to its mass.

# Law: E = m * c**2, where
## E is total (or relativistic) energy of body/system,
## m is mass (relativistic mass, not only rest mass),
## c is speed of light.

relativistic_energy = Symbol("energy", units.energy)
relativistic_mass = clone_symbol(symbols.basic.mass, "relativistic_mass")

law = Eq(relativistic_energy, relativistic_mass * units.speed_of_light**2)


@validate_input(relativistic_mass_=relativistic_mass)
@validate_output(relativistic_energy)
def calculate_rest_energy(relativistic_mass_: Quantity) -> Quantity:
    result_expr = solve(law, relativistic_energy, dict=True)[0][relativistic_energy]
    energy_applied = result_expr.subs({relativistic_mass: relativistic_mass_})
    return Quantity(energy_applied)
