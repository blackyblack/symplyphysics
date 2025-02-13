"""
Transient attenuation of separate loop coupler with cascade connection
======================================================================

The multistage connection of the couplers allows you to expand the working band and
realize less transient attenuation. Transient attenuation is a part of the power passing
into the auxiliary channel of the coupler.

..
    TODO: find link
"""

from sympy import Eq, solve, sin, asin, log
from symplyphysics import (
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

coupler_attenuation = clone_as_symbol(symbols.attenuation, subscript="0")
"""
Transient :symbols:`attenuation` of one coupler from the cascade, measured in decibels.
"""

cascade_attenuation = symbols.attenuation
"""
Transient :symbols:`attenuation` of the cascade, measured in decibels.
"""

coupler_count = symbols.positive_number
"""
Number of couplers in the cascade. See :symbols:`positive_number`.
"""

law = Eq(coupler_attenuation, 20 * log(sin(coupler_count * asin(10**(cascade_attenuation / 20))), 10))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(attenuation_of_cascade_=cascade_attenuation,
    number_of_couplers_=coupler_count)
@validate_output(coupler_attenuation)
def calculate_attenuation_of_coupler(attenuation_of_cascade_: float,
    number_of_couplers_: float) -> float:
    result_expr = solve(law, coupler_attenuation, dict=True)[0][coupler_attenuation]
    result_expr = result_expr.subs({
        cascade_attenuation: attenuation_of_cascade_,
        coupler_count: number_of_couplers_,
    })
    return convert_to_float(result_expr)
