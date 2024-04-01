from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, global_index, SumIndexed, SymbolIndexed)

# Description
## If dipoles (resistor, capacitor or coil) are connected in parallel, total admittance is a sum of admittance of each dipole.
## Law: Y_parallel = sum(Y[i]), where
## Y_parallel is total admittance,
## Y[i] is admittance of i-th dipole.

parallel_admittance = Symbol("parallel_admittance", units.conductance)
admittance = SymbolIndexed("admittance", units.conductance)
law = Eq(parallel_admittance, SumIndexed(admittance[global_index], global_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(admittances_=admittance)
@validate_output(units.conductance)
def calculate_parallel_admittance(admittances_: list[Quantity]) -> Quantity:
    local_index = Idx("local_index", (1, len(admittances_)))
    admittances_law = law.subs(global_index, local_index)
    admittances_law = admittances_law.doit()
    solved = solve(admittances_law, parallel_admittance, dict=True)[0][parallel_admittance]
    for i, v in enumerate(admittances_):
        solved = solved.subs(admittance[i + 1], v)
    return Quantity(solved)
