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
    units,
    Quantity,
    Symbol,
    dimensionless,
    validate_input,
    validate_output,
    convert_to_float,
)

relative_refractive_index = Symbol("relative_refractive_index", dimensionless)
"""
Relative refractive index between two media.

Symbol:
    :code:`n`
"""

incident_wave_speed = Symbol("incident_wave_speed", units.velocity)
r"""
Speed of the incident wave.

Symbol:
    :code:`v_incident`

Latex:
    :math:`v_\text{incident}`
"""

refracted_wave_speed = Symbol("refracted_wave_speed", units.velocity)
r"""
Speed of the refracted wave.

Symbol:
    :code:`v_refracted`

Latex:
    :math:`v_\text{refracted}`
"""

definition = Eq(relative_refractive_index, incident_wave_speed / refracted_wave_speed)
r"""
:code:`n = v_incident / v_refracted`

Latex:
    .. math::
        n = \frac{v_\text{incident}}{v_\text{refracted}}
"""


@validate_input(incident_wave_speed_=incident_wave_speed, refracted_wave_speed_=refracted_wave_speed)
@validate_output(relative_refractive_index)
def calculate_refractive_index(incident_wave_speed_: Quantity, refracted_wave_speed_: Quantity) -> float:
    result_index_expr = solve(definition, relative_refractive_index, dict=True)[0][relative_refractive_index]
    result_expr = result_index_expr.subs({
        incident_wave_speed: incident_wave_speed_,
        refracted_wave_speed: refracted_wave_speed_
    })
    result = Quantity(result_expr)
    return convert_to_float(result)
