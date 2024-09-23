"""
Mass fraction of mixture component
==================================

Mass fraction is the ratio of the mass of a mixture component to the total mass of the mixture.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    clone_as_symbol,
    symbols,
)
from symplyphysics.core.symbols.fraction import Fraction

mass_fraction = clone_as_symbol(symbols.mass_fraction, display_symbol="w[i]", display_latex="w_i")
"""
:symbols:`mass_fraction` of the mixture component.
"""

mass_of_component = clone_as_symbol(symbols.mass, display_symbol="m[i]", display_latex="m_i")
"""
:symbols:`mass` of the mixture component.
"""

mass_of_mixture = symbols.mass
"""
Total :symbols:`mass` of the mixture.
"""

definition = Eq(mass_fraction, mass_of_component / mass_of_mixture)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_of_component_=mass_of_component, mass_of_mixture_=mass_of_mixture)
@validate_output(mass_fraction)
def calculate_mass_fraction(mass_of_component_: Quantity, mass_of_mixture_: Quantity) -> Fraction:
    result_mass_fraction_expr = solve(definition, mass_fraction, dict=True)[0][mass_fraction]
    result_expr = result_mass_fraction_expr.subs({
        mass_of_component: mass_of_component_,
        mass_of_mixture: mass_of_mixture_
    })
    result = Quantity(result_expr)
    return Fraction(convert_to_float(result))
