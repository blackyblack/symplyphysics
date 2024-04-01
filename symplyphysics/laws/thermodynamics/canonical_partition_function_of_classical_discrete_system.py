from typing import Sequence
from sympy import Eq, Idx, solve
from symplyphysics import (
    dimensionless,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    SymbolIndexed,
    SumIndexed,
    global_index,
)

# Description
## Let us assume a canonical ensemble, i.e. a thermodynamically large system that is in thermal contact
## with the environment, with a temperature T and whose volume and number of constituent particles remain
## constant. For a classical discrete system the partition function is the sum of the Boltzmann factors
## of all the possible energy states:

# Law: Z = sum_i(f_i)
## Z - partition function
## f_i - Boltzmann factor of energy state i
## sum_i - summation over all energy states

partition_function = Symbol("partition_function", dimensionless)
boltzmann_factor = SymbolIndexed("boltzmann_factor", dimensionless)
law = Eq(partition_function, SumIndexed(boltzmann_factor[global_index], global_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(boltzmann_factors_=boltzmann_factor)
@validate_output(partition_function)
def calculate_partition_function(boltzmann_factors_: Sequence[float]) -> float:
    local_index = Idx("index_local", (1, len(boltzmann_factors_)))
    partition_function_law = law.subs(global_index, local_index)
    partition_function_law = partition_function_law.doit()
    solved = solve(partition_function_law, partition_function, dict=True)[0][partition_function]
    for i, v in enumerate(boltzmann_factors_):
        solved = solved.subs(boltzmann_factor[i + 1], v)
    return float(solved)
