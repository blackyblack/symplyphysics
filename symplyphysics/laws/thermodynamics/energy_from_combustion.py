from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## The amount of heat released during complete combustion of a mass
## (for solid and liquid substances) or volumetric (for gaseous) unit of matter
## Q = k_q * m
## Where:
## k_q - the heat of combustion attributed to a unit of mass or volume of fuel is called specific heat of combustion
## m - mass of matter
##
## NOTICE: This law is similar to the law of the law of energy to melt.
## Mathematically, they are similar, but different in terms of physical meaning.
## In the law of melting energy, the specific heat of melting is the amount of heat needed to melt
## a solid material weighing 1 kg. In the law of energy released during combustion,
## the specific heat of combustion is the amount of heat released during
## the complete combustion of a substance weighing 1 kg.

amount_energy = Symbol("amount_energy", units.energy)
specific_heat_combustion = Symbol("specific_heat_combustion", units.energy / units.mass)
mass_of_matter = Symbol("mass_of_matter", units.mass)

law = Eq(amount_energy, specific_heat_combustion * mass_of_matter)


def print_law() -> str:
    return print_expression(law)


@validate_input(specific_heat_combustion_=specific_heat_combustion,
    mass_of_matter_=mass_of_matter)
@validate_output(amount_energy)
def calculate_amount_energy(specific_heat_combustion_: Quantity, mass_of_matter_: Quantity) -> Quantity:
    result_amount_energy_expr = solve(law, amount_energy, dict=True)[0][amount_energy]
    result_expr = result_amount_energy_expr.subs({
        specific_heat_combustion: specific_heat_combustion_,
        mass_of_matter: mass_of_matter_
    })
    return Quantity(result_expr)
