"""
Luminosity of Sun in past from luminosity in present
====================================================

The luminosity of the Sun in the past can be calculated from the luminosity of the Sun in the present.

**Notes:**

#. The :ref:`formula for future luminosity <Luminosity of Sun in future from luminosity in present>`
   is not applicable here due to non-linearity.

..
    TODO find link
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

past_luminosity = clone_as_symbol(symbols.luminosity, subscript="0")
"""
:symbols:`luminosity` of the Sun at some point in the past.
"""

present_luminosity = symbols.luminosity
"""
:symbols:`luminosity` of the Sun in the present.
"""

time = symbols.time
"""
:symbols:`time` between the past even and today.
"""

one_billion_years = Quantity(1e9 * units.common_year,
    display_symbol="1 Gyr",
    display_latex="1 \\, \\text{Gyr}")
"""
A quantity equal to one billion years.
"""

law = Eq(past_luminosity, present_luminosity / (1 + 0.4 * (1 - ((time / one_billion_years) / 4.6))))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(luminosity_present_=present_luminosity, time_=time)
@validate_output(past_luminosity)
def calculate_luminosity_past(luminosity_present_: Quantity, time_: Quantity) -> Quantity:
    result_expr = solve(law, past_luminosity, dict=True)[0][past_luminosity]
    result_expr = result_expr.subs({
        present_luminosity: luminosity_present_,
        time: time_,
    })
    return Quantity(result_expr)
