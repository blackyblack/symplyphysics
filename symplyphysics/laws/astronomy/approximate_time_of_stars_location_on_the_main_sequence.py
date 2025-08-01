"""
Approximate lifetime of stars located on the main sequence
==========================================================

The main sequence is a stage of stellar evolution following the protostar stage. At the beginning
of the main sequence stage, the age of the star is considered to be zero. It is possible to
calculate the time spent on the main sequence for a star by knowing such a time for the Sun, as
well as the masses and luminosities of the Sun and the star.

**Notation:**

#. :quantity_notation:`sun_luminosity`.
#. :quantity_notation:`solar_mass`.

**Links:**

#. `Wikipedia, second formula <https://en.wikipedia.org/wiki/Main_sequence#Lifetime>`__.
"""

from sympy import Eq, solve
from symplyphysics import symbols, units, Quantity, validate_input, validate_output, quantities

star_lifetime = symbols.time
"""
Lifetime (:symbols:`time`) of the star on the main sequence.
"""

star_mass = symbols.mass
"""
:symbols:`mass` of the star.
"""

star_luminosity = symbols.luminosity
"""
:symbols:`luminosity` of the star.
"""

sun_lifetime = Quantity(1e10 * units.common_year, display_symbol="t_Sun", display_latex="t_\\odot")
"""
Lifetime (:symbols:`time`) of the Sun on the main sequence.
"""

law = Eq(
    star_lifetime,
    sun_lifetime * (star_mass / quantities.solar_mass) *
    (quantities.sun_luminosity / star_luminosity))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_of_star_=star_mass, luminosity_of_star_=star_luminosity)
@validate_output(star_lifetime)
def calculate_lifetime(mass_of_star_: Quantity, luminosity_of_star_: float) -> Quantity:
    result_expr = solve(law, star_lifetime, dict=True)[0][star_lifetime]
    result_expr = result_expr.subs({
        star_mass: mass_of_star_,
        star_luminosity: luminosity_of_star_,
    })
    return Quantity(result_expr)
