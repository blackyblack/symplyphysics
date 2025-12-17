"""
Kepler's constant via attracting body mass
==========================================

*Kepler's constant* is a physical quantity that only depends on the gravitational constant and the mass of
the orbited body, such as the Sun. It is constant in the sense that all planets that orbit the same body
approximately have the same value of the Kepler's constant.

**Notation:**

#. :quantity_notation:`gravitational_constant`.
"""

from sympy import Eq, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)

kepler_constant = symbols.kepler_constant
"""
:symbols:`kepler_constant` of the system.
"""

attracting_mass = clone_as_symbol(symbols.mass, display_symbol="M", display_latex="M")
"""
The :symbols:`mass` of the attracting body.
"""

law = Eq(kepler_constant, quantities.gravitational_constant * attracting_mass / (4 * pi**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(attracting_body_mass_=attracting_mass)
@validate_output(kepler_constant)
def calculate_keplers_constant(attracting_body_mass_: Quantity) -> Quantity:
    result = law.rhs.subs({
        attracting_mass: attracting_body_mass_,
    })
    return Quantity(result)
