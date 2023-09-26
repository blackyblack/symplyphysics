from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## If dipoles (resistor, capacitor or coil) are connected in parallel, total admittance is a sum of admittance of each dipole.
## Law: Y_parallel = sum(Y[i]), where
## Y_parallel is total admittance,
## Y[i] is admittance of i-th dipole.

admittances = Symbol("admittances", units.conductance)
parallel_admittance = Symbol("parallel_admittance", units.conductance)
law = Eq(parallel_admittance, SumArray(admittances), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(admittances_=admittances)
@validate_output(units.conductance)
def calculate_parallel_admittance(admittances_: list[Quantity]) -> Quantity:
    admittance_symbols = tuple_of_symbols("admittance", units.conductance, len(admittances_))
    admittances_law = law.subs(admittances, admittance_symbols).doit()
    solved = solve(admittances_law, parallel_admittance, dict=True)[0][parallel_admittance]
    for (from_, to_) in zip(admittance_symbols, admittances_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
