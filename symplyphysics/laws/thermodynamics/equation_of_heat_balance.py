from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, print_expression, Symbol, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## The first law of thermodynamics for a closed and adiabatically isolated
## system is called the heat balance equation: in a closed system of bodies,
## the algebraic sum of the amounts of heat given and received by all bodies
## involved in heat exchange is zero:
## sum(Q) = 0
## Where Q is amount of heat

amounts_energy = Symbol("amounts_energy", units.energy)
law = Eq(SumArray(amounts_energy), 0, evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(amounts_energy_=amounts_energy)
@validate_output(units.energy)
def calculate_amount_energy(amounts_energy_: list[Quantity]) -> Quantity:
    amounts_energy_symbols = tuple_of_symbols("amount_energy", units.energy, len(amounts_energy_) + 1)
    unknown_amount_energy = amounts_energy_symbols[len(amounts_energy_)]
    voltages_law = law.subs(amounts_energy, amounts_energy_symbols).doit()
    solved = solve(voltages_law, unknown_amount_energy, dict=True)[0][unknown_amount_energy]
    for (from_, to_) in zip(amounts_energy_symbols, amounts_energy_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
