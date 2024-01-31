from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, print_expression, Symbol, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols
from typing import Sequence

# Description
## The pressure of a mixture of gases is equal to the sum of their partial pressures.
## Partial pressure of a gas is the pressure that the gas would exert on the walls of a vessel while alone in it.

## Law: p = sum(p_partial)
## Where:
## p_partial is partial pressure of a gas

total_pressure = Symbol("total_pressure", units.pressure)
partial_pressures = Symbol("partial_pressures", units.pressure)

law = Eq(total_pressure, SumArray(partial_pressures), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(total_pressure_=total_pressure)
@validate_output(units.pressure)
def calculate_total_pressure(total_pressure_: Sequence[Quantity]) -> Quantity:
    partial_pressures_symbols = tuple_of_symbols("partial_pressure", units.pressure,
        len(total_pressure_))
    partial_pressures_law = law.subs(partial_pressures, partial_pressures_symbols).doit()
    solved = solve(partial_pressures_law, total_pressure, dict=True)[0][total_pressure]
    for (from_, to_) in zip(partial_pressures_symbols, total_pressure_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
