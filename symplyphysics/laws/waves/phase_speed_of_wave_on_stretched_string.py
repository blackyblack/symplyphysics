"""
Phase speed of wave on stretched string
=======================================

The phase speed of a wave on a stretched ideal string is set by properties of the string
and not by properties of the wave such as frequency and amplitude.
"""

from sympy import Eq, sqrt
from symplyphysics import (
    clone_as_symbol,
    symbols,
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
)

phase_speed = Symbol("phase_speed", units.velocity)
"""
Phase speed of the wave.

Symbol:
    :code:`v`
"""

tension = clone_as_symbol(symbols.force, display_symbol="tau", display_latex="\\tau")
"""
Tension in the string. See :symbols:`force`.
"""

linear_density = Symbol("linear_density", units.mass / units.length)
r"""
Linear density of the string, i.e. its mass per units length.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

law = Eq(phase_speed, sqrt(tension / linear_density))
r"""
:code:`v = sqrt(tau / mu)`

Latex:
    .. math::
        v = \sqrt{\frac{\tau}{\mu}}
"""

# TODO: derive from Newton's second law


@validate_input(
    string_tension_=tension,
    string_linear_density_=linear_density,
)
@validate_output(phase_speed)
def calculate_wave_speed(
    string_tension_: Quantity,
    string_linear_density_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        tension: string_tension_,
        linear_density: string_linear_density_,
    })
    return Quantity(result)
