from typing import Sequence
from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, print_expression, Symbol, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## The total number of moles in the mixture is equal to the sum of the number of moles in each of the components
## sum(nu_i) = nu
## Where:
## nu_i - moles count of mixture component
## nu - moles count of mixture
##
# Conditions:
## - Mixture in a closed impenetrable volume, that is, molecules/atoms cannot leave it
##   and they are always inside;
## - Mass is not transformed to energy, for example due to annihilation.

moles_count_of_components = Symbol("moles_count_of_components", units.amount_of_substance)
moles_count_of_mixture = Symbol("moles_count_of_mixture", units.amount_of_substance)
law = Eq(moles_count_of_mixture, SumArray(moles_count_of_components), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(moles_count_of_mixture_=moles_count_of_mixture)
@validate_output(units.amount_of_substance)
def calculate_moles_count_of_mixture(moles_count_of_mixture_: Sequence[Quantity]) -> Quantity:
    moles_count_of_component_symbols = tuple_of_symbols("moles_count_of_component", units.amount_of_substance, len(moles_count_of_mixture_))
    moles_counts_of_components_law = law.subs(moles_count_of_components, moles_count_of_component_symbols).doit()
    solved = solve(moles_counts_of_components_law, moles_count_of_mixture, dict=True)[0][moles_count_of_mixture]
    for (from_, to_) in zip(moles_count_of_component_symbols, moles_count_of_mixture_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
