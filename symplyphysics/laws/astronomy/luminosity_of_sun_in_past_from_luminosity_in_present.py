"""
Luminocity of Sun in past from luminocity in present
====================================================

The luminosity of the Sun in the past can be calculated from the luminosity of the Sun in the present.

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

past_luminosity = clone_as_symbol(symbols.luminocity, subscript="0")
"""
:symbols:`luminocity` of the Sun at some point in the past.
"""

present_luminocity = symbols.luminocity
"""
:symbols:`luminocity` of the Sun in the present.
"""

time = symbols.time
"""
:symbols:`time` between the past even and today.
"""

one_billion_years = Quantity(1e9 * units.common_year, display_symbol="1 Gyr", display_latex="1 \\, \\text{Gyr}")
"""
A quantity equal to one billion years.
"""

law = Eq(past_luminosity, present_luminocity / (1 + 0.4 * (1 - ((time / one_billion_years) / 4.6))))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(luminosity_present_=present_luminocity, time_=time)
@validate_output(past_luminosity)
def calculate_luminosity_past(luminosity_present_: Quantity, time_: Quantity) -> float:
    result_expr = solve(law, past_luminosity, dict=True)[0][past_luminosity]
    result_expr = result_expr.subs({
        present_luminocity: luminosity_present_,
        time: time_,
    })
    return Quantity(result_expr)
