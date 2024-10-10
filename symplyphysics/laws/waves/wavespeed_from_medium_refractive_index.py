"""
Wave speed from medium
======================

Speed of electromagnetic wave propagation depends on the refractive index of the medium.

**Notation:**

#. :quantity_notation:`speed_of_light`.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, quantities, symbols

wave_speed = symbols.phase_speed
"""
Speed of wave propagation in medium. See :symbols:`phase_speed`.
"""

refractive_index = symbols.relative_refractive_index
"""
:symbols:`relative_refractive_index` of the medium.
"""

law = Eq(wave_speed, quantities.speed_of_light / refractive_index)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(refraction_factor_=refractive_index)
@validate_output(wave_speed)
def calculate_wavespeed(refraction_factor_: float) -> Quantity:
    result_expr = solve(law, wave_speed, dict=True)[0][wave_speed]
    wavespeed_applied = result_expr.subs(refractive_index, refraction_factor_)
    return Quantity(wavespeed_applied)
