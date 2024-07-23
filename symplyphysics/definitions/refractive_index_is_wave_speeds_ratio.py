"""
Relative refractive index is ratio of wave speeds
=================================================

If a wave is moving from one medium to another, it refracts due to the difference
of propagation speeds in the two  media. Relative refractive index describes how much
slower the wave propagates in the refracting medium relative to the incident medium.

**Conditions:**

#. The media are isotropic and transparent.
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

refractive_index = Symbol("refractive_index", dimensionless)
"""
Relative refractive index between two media.

Symbol:
    :code:`n`
"""

outer_speed = Symbol("outer_speed", units.velocity)
r"""
Speed of the incident wave.

Symbol:
    :code:`v_incident`

Latex:
    :math:`v_\text{incident}`
"""

refracting_speed = Symbol("refracting_speed", units.velocity)
r"""
Speed of the refracted wave.

Symbol:
    :code:`v_refracted`

Latex:
    :math:`v_\text{refracted}`
"""

definition = Eq(refractive_index, outer_speed / refracting_speed)
r"""
:code:`n = v_incident / v_refracted`

Latex:
    .. math::
        n = \frac{v_\text{incident}}{v_\text{refracted}}
"""


@validate_input(outer_speed_=outer_speed, refracting_speed_=refracting_speed)
@validate_output(refractive_index)
def calculate_refractive_index(outer_speed_: Quantity, refracting_speed_: Quantity) -> float:
    result_index_expr = solve(definition, refractive_index, dict=True)[0][refractive_index]
    result_expr = result_index_expr.subs({
        outer_speed: outer_speed_,
        refracting_speed: refracting_speed_
    })
    result = Quantity(result_expr)
    return convert_to_float(result)
