from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## The heat of melting is the amount of heat that must be brought to a solid
## crystalline substance at constant pressure in order to completely transfer it
## to a liquid state.
## Q = k_lambda * m
## Where:
## Q = energy, required to melt solid matter
## k_lambda - specific heat of melting - the amount of heat that must be communicated
## to one unit of the mass of a crystalline substance in an equilibrium
## isobaric-isothermal process in order to transfer it from a solid (crystalline)
## state to a liquid one
##
## NOTICE: This law is similar to the law of energy released during
## the combustion of matter. Mathematically, they are similar,
## but different in terms of physical meaning. In the law of melting energy,
## the specific heat of melting is the amount of heat needed to melt
## a solid material weighing 1 kg. In the law of energy released during combustion,
## the specific heat of combustion is the amount of heat released during
## the complete combustion of a substance weighing 1 kg.

amount_energy = Symbol("amount_energy", units.energy)
specific_heat_melting = Symbol("specific_heat_melting", units.energy / units.mass)
mass_of_matter = Symbol("mass_of_matter", units.mass)

law = Eq(amount_energy, specific_heat_melting * mass_of_matter)


def print_law() -> str:
    return print_expression(law)


@validate_input(specific_heat_melting_=specific_heat_melting,
                mass_of_matter_=mass_of_matter)
@validate_output(amount_energy)
def calculate_amount_energy(specific_heat_melting_: Quantity, mass_of_matter_: Quantity) -> Quantity:
    result_amount_energy_expr = solve(law, amount_energy, dict=True)[0][amount_energy]
    result_expr = result_amount_energy_expr.subs({
        specific_heat_melting: specific_heat_melting_,
        mass_of_matter: mass_of_matter_
    })
    return Quantity(result_expr)
