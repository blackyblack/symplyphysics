"""
Average power of sinusoidal wave on stretched string
====================================================

The average power of a wave of any type is proportional to the square of its amplitude and to the
square of its angular frequency.

**Conditions:**

#. The wave is sinusoidal.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

wave_average_power = symbols.power
"""
Average :symbols:`power`, or rate of energy transfer, of the wave.
"""

string_linear_density = symbols.linear_density
"""
:symbols:`linear_density` of the string.
"""

phase_speed = symbols.phase_speed
"""
:symbols:`phase_speed` of the wave.
"""

wave_angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the wave.
"""

wave_amplitude = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="u_max",
    display_latex="u_\\text{max}",
)
"""
Amplitude of the wave. See :symbols:`euclidean_distance`.
"""

law = Eq(
    wave_average_power,
    string_linear_density * phase_speed * wave_angular_frequency**2 * wave_amplitude**2 / 2,
)
r"""
:code:`P = 1/2 * mu * v * w^2 * u_max^2`

Latex:
    .. math::
        P = \frac{1}{2} \mu v \omega^2 u_\text{max}^2
"""


@validate_input(
    string_linear_density_=string_linear_density,
    phase_speed_=phase_speed,
    wave_angular_frequency_=wave_angular_frequency,
    wave_amplitude_=wave_amplitude,
)
@validate_output(wave_average_power)
def calculate_average_power(
    string_linear_density_: Quantity,
    phase_speed_: Quantity,
    wave_angular_frequency_: Quantity,
    wave_amplitude_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        string_linear_density: string_linear_density_,
        phase_speed: phase_speed_,
        wave_angular_frequency: wave_angular_frequency_,
        wave_amplitude: wave_amplitude_,
    })
    return Quantity(result)
