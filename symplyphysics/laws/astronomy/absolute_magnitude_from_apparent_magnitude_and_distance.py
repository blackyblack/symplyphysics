"""
Absolute magnitude of stars from apparent magnitude and distance
================================================================

Absolute magnitude can be calculated using apparent magnitude and distance to the object.
"""

from sympy import Eq, solve, log
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

absolute_magnitude = SymbolNew("M", dimensionless)
"""
Absolute magnitude for stars is defined as the apparent magnitude of an object if it were located at a distance
of 10 parsecs (2.063e+6 astronomical units) from the observer and would not experience either interstellar or atmospheric absorption.
"""

apparent_magnitude = SymbolNew("m", dimensionless)
"""
The apparent magnitude is a measure of the brightness of a celestial body (more precisely, the illumination created
by this body) from the observer's point of view. The brighter the object, the smaller its magnitude.
"""

distance = SymbolNew("d", units.length)
"""
Distance to the object.
"""

distance_constant = Quantity(2.063e+6 * units.astronomical_unit,
    display_symbol="d0",
    display_latex="d_0")
"""
Constant equal to 2.063e+6 astronomical units.
"""

law = Eq(absolute_magnitude, apparent_magnitude - 5 * log(distance / distance_constant, 10))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(apparent_magnitude_=apparent_magnitude, distance_=distance)
@validate_output(absolute_magnitude)
def calculate_absolute_magnitude(apparent_magnitude_: float, distance_: Quantity) -> float:
    result_expr = solve(law, absolute_magnitude, dict=True)[0][absolute_magnitude]
    result_expr = result_expr.subs({
        apparent_magnitude: apparent_magnitude_,
        distance: distance_,
    })
    return convert_to_float(result_expr)
