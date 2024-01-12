from typing import Sequence
from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, print_expression, Symbol, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## The mass of a mixture of liquids (gases) is equal to the sum of the masses of the components of the mixture
## sum(m_i) = m
## Where:
## m_i - mass of mixture component
## m - mass of mixture

masses_of_components = Symbol("masses_of_components", units.mass)
mass_of_mixture = Symbol("mass_of_mixture", units.mass)
law = Eq(mass_of_mixture, SumArray(masses_of_components), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_of_mixture_=mass_of_mixture)
@validate_output(units.mass)
def calculate_mass_of_mixture(mass_of_mixture_: Sequence[Quantity]) -> Quantity:
    mass_of_component_symbols = tuple_of_symbols("mass_of_component", units.force, len(mass_of_mixture_))
    masses_of_components_law = law.subs(masses_of_components, mass_of_component_symbols).doit()
    solved = solve(masses_of_components_law, mass_of_mixture, dict=True)[0][mass_of_mixture]
    for (from_, to_) in zip(mass_of_component_symbols, mass_of_mixture_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
