from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols)

# Description
## The heat of vaporization of a substance is the amount of heat that must be communicated
# to a substance at constant pressure and temperature in order to transfer
# it from a liquid state to a gaseous state (into steam).
## Q = k_v * m
## Where:
## Q - the energy required to convert a substance from a liquid to a gaseous state
## k_v - the specific heat of vaporization is the necessary energy that must be communicated
## to a substance weighing 1 kg in order to convert it to a gaseous state
##
## NOTICE: This law is similar to the law of energy released during
## the combustion of matter. Mathematically, they are similar,
## but different in terms of physical meaning. In the law of melting energy,
## the specific heat of melting is the amount of heat needed to melt
## a solid material weighing 1 kg. In the law of energy released during combustion,
## the specific heat of combustion is the amount of heat released during
## the complete combustion of a substance weighing 1 kg.
## In the law of energy for vaporization, the specific heat of vaporization is the amount of energy
## must be expended to evaporate one kilogram of a substance taken at boiling point.

amount_energy = Symbol("amount_energy", units.energy)
mass = symbols.basic.mass
specific_heat_vaporization = Symbol("specific_heat_vaporization", units.energy / units.mass)

law = Eq(amount_energy, specific_heat_vaporization * mass)


def print_law() -> str:
    return print_expression(law)


@validate_input(specific_heat_vaporization_=specific_heat_vaporization,
    mass_of_matter_=mass)
@validate_output(amount_energy)
def calculate_amount_energy(specific_heat_vaporization_: Quantity,
    mass_of_matter_: Quantity) -> Quantity:
    result_amount_energy_expr = solve(law, amount_energy, dict=True)[0][amount_energy]
    result_expr = result_amount_energy_expr.subs({
        specific_heat_vaporization: specific_heat_vaporization_,
        mass: mass_of_matter_
    })
    return Quantity(result_expr)
