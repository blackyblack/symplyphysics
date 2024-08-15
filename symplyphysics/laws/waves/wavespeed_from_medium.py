"""
Wave speed from medium
======================

Speed of electromagnetic wave propagation depends on the refractive index of the medium.
"""

from sympy.physics.units import speed_of_light
from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, dimensionless, validate_input,
    validate_output)

wave_speed = Symbol("wave_speed", units.velocity)
"""
Speed of wave propagation in medium.

Symbol:
    :code:`v`
"""

refractive_index = Symbol("refractive_index", dimensionless)
"""
Refractive index of the medium relative to vacuum.

Symbol:
    :code:`n`
"""

law = Eq(wave_speed, speed_of_light / refractive_index)
r"""
:code:`v = c / n`

Latex:
    .. math::
        v = \frac{c}{n}
"""


@validate_input(refraction_factor_=refractive_index)
@validate_output(wave_speed)
def calculate_wavespeed(refraction_factor_: float) -> Quantity:
    result_expr = solve(law, wave_speed, dict=True)[0][wave_speed]
    wavespeed_applied = result_expr.subs(refractive_index, refraction_factor_)
    return Quantity(wavespeed_applied)
