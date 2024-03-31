from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_indexed import SumIndexed
from symplyphysics.core.symbols.symbols import SymbolIndexed


# Description
## If dipoles (resistor, capacitor or coil) are connected in parallel, total admittance is a sum of admittance of each dipole.
## Law: Y_parallel = sum(Y[i]), where
## Y_parallel is total admittance,
## Y[i] is admittance of i-th dipole.

parallel_admittance = Symbol("parallel_admittance", units.conductance)
admittance = SymbolIndexed("admittance", units.conductance)
admittance_index = Idx("admittance_index")
law = Eq(parallel_admittance, SumIndexed(admittance[admittance_index], admittance_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(admittances_=admittance)
@validate_output(units.conductance)
def calculate_parallel_admittance(admittances_: list[Quantity]) -> Quantity:
    index_local = Idx("index_local", (1, len(admittances_)))
    admittances_law = law.subs(admittance_index, index_local)
    admittances_law = admittances_law.doit()
    solved = solve(admittances_law, parallel_admittance, dict=True)[0][parallel_admittance]
    for (from_, to_) in zip(range(1, len(admittances_) + 1), admittances_):
        solved = solved.subs(admittance[from_], to_)
    return Quantity(solved)
