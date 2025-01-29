"""
Ratio of luminosities from ratio of masses of stars
===================================================

The comparisons of masses and luminosities for most stars revealed the following relationship:
luminosity is approximately proportional to the fourth power of mass.

**Conditions:**

#. The masses must obey the inequality :math:`0.43 m_1 < m_2 < 2 m_1`, see link for more information.
   The symbols :math:`m_1, m_2` refer to those defined below.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Mass%E2%80%93luminosity_relation>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

first_mass = clone_as_symbol(symbols.mass, subscript="1")
"""
:symbols:`mass` of the first star.
"""

second_mass = clone_as_symbol(symbols.mass, subscript="2")
"""
:symbols:`mass` of the second star.
"""

first_luminosity = clone_as_symbol(symbols.luminosity, subscript="1")
"""
:symbols:`luminosity` of the first star.
"""

second_luminosity = clone_as_symbol(symbols.luminosity, subscript="2")
"""
:symbols:`luminosity` of the second star.
"""

law = Eq(second_luminosity / first_luminosity, (second_mass / first_mass)**4)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_first_=first_mass,
    mass_second_=second_mass,
    illuminance_first_=first_luminosity)
@validate_output(second_luminosity)
def calculate_illuminance_second(mass_first_: Quantity, mass_second_: Quantity,
    illuminance_first_: Quantity) -> Quantity:
    result_expr = solve(law, second_luminosity, dict=True)[0][second_luminosity]
    result_expr = result_expr.subs({
        first_mass: mass_first_,
        second_mass: mass_second_,
        first_luminosity: illuminance_first_
    })
    return Quantity(result_expr)
