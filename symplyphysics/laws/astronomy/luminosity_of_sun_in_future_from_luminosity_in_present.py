"""
Luminosity of Sun in future from luminosity in present
======================================================

The luminosity of the Sun in the future can be calculated from the luminosity of the Sun in the present.

**Notes:**

#. The :ref:`formula for luminosity in the past <Luminosity of Sun in past from luminosity in present>`
   is not applicable here due to non-linearity.

**Conditions:**

#. The formula is valid while the Sun is `on the main sequence <https://faculty.wcas.northwestern.edu/infocom/The%20Website/end.html>`__.

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

future_luminosity = clone_as_symbol(symbols.luminosity, subscript="1")
"""
:symbols:`luminosity` of the Sun at some point in the future.
"""

present_luminosity = symbols.luminosity
"""
:symbols:`luminosity` of the Sun in the present.
"""

time = symbols.time
"""
:symbols:`time` between now and the future event.
"""

one_billion_years = Quantity(1e9 * units.common_year, display_symbol="1 Gyr", display_latex="1 \\, \\text{Gyr}")
"""
A quantity equal to one billion years.
"""

law = Eq(
    future_luminosity,
    present_luminosity * ((5.59 / (time / one_billion_years)) - 1.41 + 0.26 *
    (time / one_billion_years)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(luminosity_present_=present_luminosity, time_=time)
@validate_output(future_luminosity)
def calculate_luminosity_future(luminosity_present_: Quantity, time_: Quantity) -> Quantity:
    result_expr = solve(law, future_luminosity, dict=True)[0][future_luminosity]
    result_expr = result_expr.subs({
        present_luminosity: luminosity_present_,
        time: time_,
    })
    return Quantity(result_expr)
