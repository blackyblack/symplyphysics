from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## If resistors are connected in parallel, total conductance is a sum of conductances of each resistor.
## Law: sigma_parallel = sum(sigma[i]), where
## sigma_parallel is total conductance,
## sigma[i] is conductance of i-th resistor.

conductances = Symbol("conductances", units.conductance)
parallel_conductance = Symbol("parallel_conductance", units.conductance)
law = Eq(parallel_conductance, SumArray(conductances), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(conductances_=conductances)
@validate_output(units.conductance)
def calculate_parallel_conductance(conductances_: list[Quantity]) -> Quantity:
    conductance_symbols = tuple_of_symbols("conductance", units.conductance, len(conductances_))
    currents_law = law.subs(conductances, conductance_symbols).doit()
    solved = solve(currents_law, parallel_conductance, dict=True)[0][parallel_conductance]
    for (from_, to_) in zip(conductance_symbols, conductances_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
