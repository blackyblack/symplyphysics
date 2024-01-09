from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## If resistors are connected in parallel, total resistance is the inverse sum of inverse resistance of each resistor.
## Law: R_parallel = 1/sum(1/R[i]), where
## R_parallel is total resistance,
## R[i] is resistance of i-th resistor.

inverse_resistances = Symbol("inverse_resistances", 1 / units.impedance)
parallel_resistance = Symbol("parallel_resistance", units.impedance)
law = Eq(parallel_resistance, 1 / SumArray(inverse_resistances), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(inverse_resistances_=inverse_resistances)
@validate_output(units.impedance)
def calculate_parallel_resistance(inverse_resistances_: list[Quantity]) -> Quantity:
    inverse_resistance_symbols = tuple_of_symbols("inverse_resistance", 1 / units.impedance, len(inverse_resistances_))
    resistances_law = law.subs(inverse_resistances, inverse_resistance_symbols).doit()
    solved = solve(resistances_law, parallel_resistance, dict=True)[0][parallel_resistance]
    for from_, to_ in zip(inverse_resistance_symbols, inverse_resistances_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
