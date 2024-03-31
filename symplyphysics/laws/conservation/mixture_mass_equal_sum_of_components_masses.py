from typing import Sequence
from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, print_expression, validate_input, validate_output,
    symbols, clone_symbol)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

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

masses_of_components = clone_symbol(symbols.basic.mass, "masses_of_components")
law = Eq(symbols.basic.mass, SumArray(masses_of_components), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(masses_of_components_=masses_of_components)
@validate_output(symbols.basic.mass)
def calculate_mass_of_mixture(masses_of_components_: Sequence[Quantity]) -> Quantity:
    mass_of_component_symbols = tuple_of_symbols("mass_of_component", units.mass,
        len(masses_of_components_))
    masses_of_components_law = law.subs(masses_of_components, mass_of_component_symbols).doit()
    solved = solve(masses_of_components_law, symbols.basic.mass, dict=True)[0][symbols.basic.mass]
    for (from_, to_) in zip(mass_of_component_symbols, masses_of_components_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
