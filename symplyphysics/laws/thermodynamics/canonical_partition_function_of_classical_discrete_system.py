from typing import Collection
from sympy import Eq
from symplyphysics import (
    dimensionless,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## Let us assume a canonical ensemble, i.e. a thermodynamically large system that is in thermal contact
## with the environment, with a temperature T and whose volume and number of constituent particles remain
## constant. For a classical discrete system the partition function is the sum of the Boltzmann factors
## of all the possible energy states:

# Law: Z = sum_i(f_i)
## Z - partition function
## f_i - Boltzmann factor of energy state i
## sum_i - summation over all energy states

boltzmann_factors = Symbol("boltzmann_factors", dimensionless)
partition_function = Symbol("partition_function", dimensionless)

law = Eq(partition_function, SumArray(boltzmann_factors))


def print_law() -> str:
    return print_expression(law)


@validate_input(boltzmann_factors_=boltzmann_factors)
@validate_output(partition_function)
def calculate_partition_function(boltzmann_factors_: Collection[float]) -> float:
    boltzmann_factor_symbols = tuple_of_symbols("boltzmann_factor", dimensionless, len(boltzmann_factors_))
    expr = law.rhs.subs(boltzmann_factors, boltzmann_factor_symbols).doit()
    for symbol, value in zip(boltzmann_factor_symbols, boltzmann_factors_):
        expr = expr.subs(symbol, value)
    return float(expr)
