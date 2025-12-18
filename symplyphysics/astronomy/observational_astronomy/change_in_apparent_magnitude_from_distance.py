"""
Change in apparent magnitude from distance
==========================================

The apparent magnitude is a measure of the brightness of a celestial body (more
precisely, the illumination created by this body) from the observer's point of view. The
brighter the object, the smaller its magnitude. The relationship of the stellar magnitude
scale with real physical quantities is logarithmic, since a change in brightness by the
same number of times is perceived by the eye as a change by the same amount. The
difference in the stellar magnitudes of two objects is equal to the decimal logarithm of
the ratio of their illuminances, up to a multiplier.

**Links:**

#. `Wikipedia, second formula <https://en.wikipedia.org/wiki/Apparent_magnitude#Calculations>`__.
"""

from sympy import Eq, solve, log
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    clone_as_symbol,
    symbols,
)

first_apparent_magnitude = clone_as_symbol(symbols.apparent_magnitude, subscript="1")
"""
:symbols:`apparent_magnitude` of the first object.
"""

second_apparent_magnitude = clone_as_symbol(symbols.apparent_magnitude, subscript="2")
"""
:symbols:`apparent_magnitude` of the second object.
"""

first_irradiance = clone_as_symbol(symbols.irradiance,
    display_symbol="E_e1",
    display_latex="E_{\\text{e}1}")
"""
Observed :symbols:`irradiance` of the first object.
"""

second_irradiance = clone_as_symbol(symbols.irradiance,
    display_symbol="E_e2",
    display_latex="E_{\\text{e}2}")
"""
Observed :symbols:`irradiance` of the second object.
"""

law = Eq(second_apparent_magnitude - first_apparent_magnitude,
    -2.5 * log(second_irradiance / first_irradiance, 10))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(apparent_magnitude_first_=first_apparent_magnitude,
    illuminance_first_=first_irradiance,
    illuminance_second_=second_irradiance)
@validate_output(second_apparent_magnitude)
def calculate_apparent_magnitude_second(apparent_magnitude_first_: float,
    illuminance_first_: Quantity, illuminance_second_: Quantity) -> float:
    result_expr = solve(law, second_apparent_magnitude, dict=True)[0][second_apparent_magnitude]
    result_expr = result_expr.subs({
        first_apparent_magnitude: apparent_magnitude_first_,
        first_irradiance: illuminance_first_,
        second_irradiance: illuminance_second_
    })
    return convert_to_float(result_expr)
