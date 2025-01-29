"""
Latitude from zenith angle and declination
==========================================

Knowing the zenith angles and declinations of the northern and southern stars, it is possible
to determine the latitude of the observation site.

**Notes:**

#. Northern star is any star north of the zenith with known declination.
#. Southern star is any star south of the zenith with known declination.

**Conditions:**

#. Both stars are at upper transit (culmination).

..
    TODO find link
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

latitude = symbols.latitude
"""
:symbols:`latitude` of the observation site.
"""

north_zenith_angle = clone_as_symbol(symbols.zenith_angle,
    display_symbol="theta_N",
    display_latex="\\theta_\\text{N}")
"""
:symbols:`zenith_angle` of the northern star.
"""

south_zenith_angle = clone_as_symbol(symbols.zenith_angle,
    display_symbol="theta_S",
    display_latex="\\theta_\\text{S}")
"""
:symbols:`zenith_angle` of the southern star.
"""

north_declination = clone_as_symbol(symbols.declination,
    display_symbol="delta_N",
    display_latex="\\delta_\\text{N}")
"""
:symbols:`declination` of the northern star.
"""

south_declination = clone_as_symbol(symbols.declination,
    display_symbol="delta_S",
    display_latex="\\delta_\\text{S}")
"""
:symbols:`declination` of the southern star.
"""

law = Eq(latitude,
    (south_zenith_angle - north_zenith_angle + south_declination + north_declination) / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    zenith_distance_north_=north_zenith_angle,
    south_zenith_angle_=south_zenith_angle,
    northern_declination_=north_declination,
    south_declination_=south_declination,
)
@validate_output(latitude)
def calculate_latitude(zenith_distance_north_: Quantity, south_zenith_angle_: Quantity,
    northern_declination_: Quantity, south_declination_: Quantity) -> Quantity:
    result_expr = solve(law, latitude, dict=True)[0][latitude]
    result_expr = result_expr.subs({
        north_zenith_angle: zenith_distance_north_,
        south_zenith_angle: south_zenith_angle_,
        north_declination: northern_declination_,
        south_declination: south_declination_
    })
    return Quantity(result_expr)
