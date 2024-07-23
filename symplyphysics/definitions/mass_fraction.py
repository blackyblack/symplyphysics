"""
Mass fraction of mixture component
==================================

Mass fraction is the ratio of the mass of a mixture component to the total mass of the mixture.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    clone_symbol,
    symbols,
)
from symplyphysics.core.symbols.fraction import Fraction

mass_fraction = Symbol("mass_fraction", dimensionless)
r"""
Mass fraction of the mixture component.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

mass_of_component = clone_symbol(symbols.basic.mass, "mass_of_component")
"""
:attr:`~symplyphysics.symbols.basic.mass` of the mixture component.

Symbol:
    :code:`m_i`

Latex:
    :math:`m_i`
"""

mass_of_mixture = clone_symbol(symbols.basic.mass, "mass_of_mixture")
"""
Total :attr:`~symplyphysics.symbols.basic.mass` of the mixture.

Symbol:
    :code:`m`
"""

definition = Eq(mass_fraction, mass_of_component / mass_of_mixture)
r"""
:code:`w = m_i / m`

Latex:
    .. math::
        \omega = \frac{m_i}{m}
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
