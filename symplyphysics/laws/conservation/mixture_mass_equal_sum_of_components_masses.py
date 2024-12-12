from typing import Sequence
from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, validate_input, validate_output,
    symbols, SymbolIndexed, SumIndexed, global_index)

# Description
## The mass of a mixture of liquids (gases) is equal to the sum of the masses of the components of the mixture
## sum(m_i) = m
## Where:
## m_i - mass of mixture component
## m - mass of mixture
##
# Conditions:
## - Mixture in a closed impenetrable volume, that is, molecules/atoms cannot leave it
##   and they are always inside;
## - Mass is not transformed to energy, for example due to annihilation.

# Links: Engineering LibreTexts, Composition on a Mass Basis <https://eng.libretexts.org/Bookshelves/Introductory_Engineering/Basic_Engineering_Science_-_A_Systems_Accounting_and_Modeling_Approach_(Richards)/03%3A_Conservation_of_Mass/3.04%3A_Mixture_Composition>

mass_of_mixture = symbols.mass
# TODO: clone from symbols.mass
mass_of_component = SymbolIndexed("mass_of_component", units.mass)
law = Eq(mass_of_mixture, SumIndexed(mass_of_component[global_index], global_index))


@validate_input(masses_of_components_=mass_of_component)
@validate_output(mass_of_mixture)
def calculate_mass_of_mixture(masses_of_components_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(masses_of_components_)))
    masses_of_components_law = law.subs(global_index, local_index)
    masses_of_components_law = masses_of_components_law.doit()
    solved = solve(masses_of_components_law, mass_of_mixture, dict=True)[0][mass_of_mixture]
    for i, v in enumerate(masses_of_components_):
        solved = solved.subs(mass_of_component[i + 1], v)
    return Quantity(solved)
