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

inv_resistances = Symbol("resistances", 1 / units.impedance)
parallel_resistance = Symbol("parallel_resistance", units.impedance)
law = Eq(parallel_resistance, 1 / SumArray(inv_resistances), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(resistances_=units.impedance)
@validate_output(units.impedance)
def calculate_parallel_resistance(resistances_: list[Quantity]) -> Quantity:
    resistance_symbols = tuple_of_symbols("resistance", units.impedance, len(resistances_))
    inv_resistance_symbols = tuple(1 / r for r in resistance_symbols)
    resistances_law = law.subs(inv_resistances, inv_resistance_symbols).doit()
    solved = solve(resistances_law, parallel_resistance, dict=True)[0][parallel_resistance]
    for from_, to_ in zip(resistance_symbols, resistances_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
