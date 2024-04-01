from typing import Sequence
from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, print_expression, Symbol, validate_input,
    validate_output, SymbolIndexed, SumIndexed, global_index)

# Description
## The pressure of a mixture of gases is equal to the sum of their partial pressures.
## Partial pressure of a gas is the pressure that the gas would exert on the walls of a vessel while alone in it.

## Law: p = sum(p_partial)
## Where:
## p_partial is partial pressure of a gas

total_pressure = Symbol("total_pressure", units.pressure)
partial_pressure = SymbolIndexed("partial_pressure", units.pressure)
law = Eq(total_pressure, SumIndexed(partial_pressure[global_index], global_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(partial_pressures_=partial_pressure)
@validate_output(total_pressure)
def calculate_total_pressure(partial_pressures_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(partial_pressures_)))
    partial_pressures_law = law.subs(global_index, local_index)
    partial_pressures_law = partial_pressures_law.doit()
    solved = solve(partial_pressures_law, total_pressure, dict=True)[0][total_pressure]
    for i, v in enumerate(partial_pressures_):
        solved = solved.subs(partial_pressure[i + 1], v)
    return Quantity(solved)
