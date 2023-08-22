from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## If capacitors are connected in parallel, total capacitance is a sum of capacitances of each capacitor.
## Law: C_parallel = sum(C[i]), where
## C_parallel is total capacitance,
## C[i] is capacitance of i-th capacitor.

capacitances = Symbol("capacitances", units.capacitance)
parallel_capacitance = Symbol("parallel_capacitance", units.capacitance)
law = Eq(parallel_capacitance, SumArray(capacitances), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(capacitances_=capacitances)
@validate_output(units.capacitance)
def calculate_parallel_capacitance(capacitances_: list[Quantity]) -> Quantity:
    capacitance_symbols = tuple_of_symbols("capacitance", units.capacitance, len(capacitances_))
    capacitances_law = law.subs(capacitances, capacitance_symbols).doit()
    solved = solve(capacitances_law, parallel_capacitance, dict=True)[0][parallel_capacitance]
    for (from_, to_) in zip(capacitance_symbols, capacitances_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
