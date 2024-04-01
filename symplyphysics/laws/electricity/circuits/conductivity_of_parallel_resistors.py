from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, SymbolIndexed, SumIndexed, global_index)

# Description
## If resistors are connected in parallel, total conductance is a sum of conductances of each resistor.
## Law: sigma_parallel = sum(sigma[i]), where
## sigma_parallel is total conductance,
## sigma[i] is conductance of i-th resistor.

parallel_conductance = Symbol("parallel_conductance", units.conductance)
conductance = SymbolIndexed("conductance", units.conductance)
law = Eq(parallel_conductance, SumIndexed(conductance[global_index], global_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(conductances_=conductance)
@validate_output(parallel_conductance)
def calculate_parallel_conductance(conductances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(conductances_)))
    conductances_law = law.subs(global_index, local_index)
    conductances_law = conductances_law.doit()
    solved = solve(conductances_law, parallel_conductance, dict=True)[0][parallel_conductance]
    for i, v in enumerate(conductances_):
        solved = solved.subs(conductance[i + 1], v)
    return Quantity(solved)
