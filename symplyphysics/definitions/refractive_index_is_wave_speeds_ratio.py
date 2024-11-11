"""
Relative refractive index is ratio of wave speeds
=================================================

If a wave is moving from one medium to another, it refracts due to the difference
of propagation speeds in the two  media. Relative refractive index describes how much
slower the wave propagates in the refracting medium relative to the incident medium.

**Conditions:**

#. Both media are isotropic and transparent.
#. The wave is monochromatic. Note that the speed of wave propagation depends on the
   wave frequency.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

relative_refractive_index = symbols.relative_refractive_index
"""
:symbols:`relative_refractive_index` between the two media.
"""

incident_wave_speed = clone_as_symbol(
    symbols.phase_speed,
    display_symbol="v_incident",
    display_latex="v_\\text{incident}",
)
"""
:symbols:`phase_speed` of the incident wave.
"""

refracted_wave_speed = clone_as_symbol(
    symbols.phase_speed,
    display_symbol="v_refracted",
    display_latex="v_\\text{refracted}",
)
"""
:symbols:`phase_speed` of the refracted wave.
"""

definition = Eq(relative_refractive_index, incident_wave_speed / refracted_wave_speed)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(incident_wave_speed_=incident_wave_speed,
    refracted_wave_speed_=refracted_wave_speed)
@validate_output(relative_refractive_index)
def calculate_refractive_index(incident_wave_speed_: Quantity,
    refracted_wave_speed_: Quantity) -> float:
    result_index_expr = solve(definition, relative_refractive_index,
        dict=True)[0][relative_refractive_index]
    result_expr = result_index_expr.subs({
        incident_wave_speed: incident_wave_speed_,
        refracted_wave_speed: refracted_wave_speed_
    })
    result = Quantity(result_expr)
    return convert_to_float(result)
