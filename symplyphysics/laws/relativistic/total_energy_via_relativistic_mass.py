from sympy import Eq, solve
from symplyphysics import units, Quantity, Symbol, validate_input, validate_output, symbols

# Description
## Fundamentally the energy of an object is synonimical to its mass.

# Law: E = m * c**2, where
## E is total (or relativistic) energy of body/system,
## m is mass (relativistic mass, not only rest mass),
## c is speed of light.

energy = Symbol("energy", units.energy)
mass = symbols.basic.mass

law = Eq(energy, mass * units.speed_of_light**2)


@validate_input(mass_=mass)
@validate_output(energy)
def calculate_rest_energy(mass_: Quantity) -> Quantity:
    result_expr = solve(law, energy, dict=True)[0][energy]
    energy_applied = result_expr.subs({mass: mass_})
    return Quantity(energy_applied)
