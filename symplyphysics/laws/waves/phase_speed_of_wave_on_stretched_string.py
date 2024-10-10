"""
Phase speed of wave on stretched string
=======================================

The phase speed of a wave on a stretched ideal string is set by properties of the string
and not by properties of the wave such as frequency and amplitude.
"""

from sympy import Eq, sqrt
from symplyphysics import symbols, Quantity, validate_input, validate_output

phase_speed = symbols.phase_speed
"""
:symbols:`phase_speed` of the wave.
"""

tension = symbols.tension
"""
:symbols:`tension` in the string.
"""

linear_density = symbols.linear_density
"""
:symbols:`linear_density` of the string, i.e. its mass per units length.
"""

law = Eq(phase_speed, sqrt(tension / linear_density))
"""
:laws:symbol::

:laws:latex::
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
