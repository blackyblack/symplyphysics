from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, SymbolIndexed, SumIndexed, global_index)

# Description
## If capacitors are connected in parallel, total capacitance is a sum of capacitances of each capacitor.
## Law: C_parallel = sum(C[i]), where
## C_parallel is total capacitance,
## C[i] is capacitance of i-th capacitor.

parallel_capacitance = Symbol("parallel_capacitance", units.capacitance)
capacitance = SymbolIndexed("capacitance", units.capacitance)
law = Eq(parallel_capacitance, SumIndexed(capacitance[global_index], global_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(capacitances_=capacitance)
@validate_output(parallel_capacitance)
def calculate_parallel_capacitance(capacitances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(capacitances_)))
    capacitances_law = law.subs(global_index, local_index)
    capacitances_law = capacitances_law.doit()
    solved = solve(capacitances_law, parallel_capacitance, dict=True)[0][parallel_capacitance]
    for i, v in enumerate(capacitances_):
        solved = solved.subs(capacitance[i + 1], v)
    return Quantity(solved)
