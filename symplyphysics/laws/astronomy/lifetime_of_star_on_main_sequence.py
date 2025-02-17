"""
Lifetime of star on main sequence
=================================

The main sequence is the stage of stellar evolution. This stage begins in stars after the
protostar stage. At the beginning of the main sequence stage, the age of the star is
considered to be zero. The lifetime of the star can be approximated from the mass of the
star.

**Notation:**

#. :quantity_notation:`solar_mass`.

**Notes:**

#. The indicator has the value of :math:`4.75` for stars with mass of :math:`0.7` to
   :math:`2` solar masses, and :math:`4.75 m + 2.125` for stars with mass of :math:`0.1`
   to :math:`0.7` solar masses. Here :math:`m` refers to the mass of the star.

..
    TODO find link
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    quantities,
)

lifetime = symbols.time
"""
The lifetime (:symbols:`time`) of the star.
"""

mass = symbols.mass
"""
:symbols:`mass` of the star.
"""

indicator = Symbol("n", dimensionless)
"""
A dimensionless parameter whose value depends on the mass of the star.
"""

ten_billion_years = Quantity(1e10 * units.common_year,
    display_symbol="10 Gyr",
    display_latex="10 \\, \\text{Gyr}")
"""
A quantity equal to ten billion years.
"""

law = Eq(lifetime, ten_billion_years * ((mass / quantities.solar_mass)**(1 - indicator)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_of_star_=mass, indicator_=indicator)
@validate_output(lifetime)
def calculate_lifetime(mass_of_star_: Quantity, indicator_: float) -> Quantity:
    result_expr = solve(law, lifetime, dict=True)[0][lifetime]
    result_expr = result_expr.subs({
        mass: mass_of_star_,
        indicator: indicator_,
    })
    return Quantity(result_expr)
