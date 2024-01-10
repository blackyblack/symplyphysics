from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description

amount_energy = Symbol("", units.energy)
specific_heat_melting = Symbol("", units.energy / units.mass)
mass_of_matter = Symbol("", units.mass)

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
