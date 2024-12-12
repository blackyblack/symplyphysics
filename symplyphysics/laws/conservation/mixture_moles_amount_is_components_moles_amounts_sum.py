from typing import Sequence
from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, SymbolIndexed, SumIndexed, global_index)

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

# Links: Engineering LibreTexts, Composition on a Molar Basis <https://eng.libretexts.org/Bookshelves/Introductory_Engineering/Basic_Engineering_Science_-_A_Systems_Accounting_and_Modeling_Approach_(Richards)/03%3A_Conservation_of_Mass/3.04%3A_Mixture_Composition>

moles_count_of_mixture = Symbol("moles_count_of_mixture", units.amount_of_substance)
moles_count_of_component = SymbolIndexed("moles_count_of_component", units.amount_of_substance)
law = Eq(moles_count_of_mixture, SumIndexed(moles_count_of_component[global_index], global_index))


@validate_input(moles_count_of_components_=moles_count_of_component)
@validate_output(moles_count_of_mixture)
def calculate_moles_count_of_mixture(moles_count_of_components_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(moles_count_of_components_)))
    moles_counts_of_mixture_law = law.subs(global_index, local_index)
    moles_counts_of_mixture_law = moles_counts_of_mixture_law.doit()
    solved = solve(moles_counts_of_mixture_law, moles_count_of_mixture,
        dict=True)[0][moles_count_of_mixture]
    for i, v in enumerate(moles_count_of_components_):
        solved = solved.subs(moles_count_of_component[i + 1], v)
    return Quantity(solved)
