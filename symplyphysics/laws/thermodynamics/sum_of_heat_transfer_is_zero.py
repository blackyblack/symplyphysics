from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, print_expression, validate_input, validate_output,
    SymbolIndexed, SumIndexed, global_index)

# Description
## The first law of thermodynamics for a closed and adiabatically isolated
## system is called the heat balance equation: in a closed system of bodies,
## the algebraic sum of the amounts of heat given and received by all bodies
## involved in heat exchange is zero:
## sum(Q) = 0
## Where Q is amount of heat

amount_energy = SymbolIndexed("amount_energy", units.energy)
law = Eq(SumIndexed(amount_energy[global_index], global_index), 0)


def print_law() -> str:
    return print_expression(law)


@validate_input(amounts_energy_=amount_energy)
@validate_output(units.energy)
def calculate_amount_energy(amounts_energy_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(amounts_energy_) + 1))
    amounts_energy_law = law.subs(global_index, local_index)
    amounts_energy_law = amounts_energy_law.doit()
    unknown_amount_energy = amount_energy[len(amounts_energy_) + 1]
    solved = solve(amounts_energy_law, unknown_amount_energy, dict=True)[0][unknown_amount_energy]
    for i, v in enumerate(amounts_energy_):
        solved = solved.subs(amount_energy[i + 1], v)
    return Quantity(solved)
