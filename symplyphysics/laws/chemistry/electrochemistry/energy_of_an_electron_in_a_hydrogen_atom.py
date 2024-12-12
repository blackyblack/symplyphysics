from sympy.physics.units import elementary_charge
from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

# Description
## There is an expression for the total energy of the hydrogen atom according to Bohr's theory.

# Law: E = k * e^2 / (2 * r), where
## E - total energy of electron,
## k - Coulomb's constant,
## e - charge of electron,
## r - radius of the Bohr orbit of electron.

# Links: derivable from first formula and formula for kinetic energy <https://en.wikipedia.org/wiki/Classical_electron_radius#>

energy_of_electron = Symbol("energy_of_electron", units.energy)

radius_of_electron = Symbol("radius_of_electron", units.length)

law = Eq(energy_of_electron,
    units.coulomb_constant * elementary_charge**2 / (2 * radius_of_electron))


@validate_input(radius_of_electron_=radius_of_electron)
@validate_output(energy_of_electron)
def calculate_energy_of_electron(radius_of_electron_: Quantity) -> Quantity:
    result_expr = solve(law, energy_of_electron, dict=True)[0][energy_of_electron]
    result = result_expr.subs(radius_of_electron, radius_of_electron_)
    return Quantity(result)
