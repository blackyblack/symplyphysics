"""
Average power of sinusoidal wave on stretched string
====================================================

The average power of a wave of any type is proportional to the square of its amplitude and to the
square of its angular frequency.
"""

from sympy import Eq
from symplyphysics import (
    units,
    angle_type,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
)

wave_average_power = Symbol("wave_average_power", units.power)
"""
Average power, or rate of energy transfer, of the wave.

Symbol:
    :code:`P`
"""

string_linear_density = Symbol("string_linear_density", units.mass / units.length)
r"""
Linear density of the string, i.e. mass per unit length.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

wave_phase_velocity = Symbol("wave_phase_velocity", units.velocity)
"""
Phase velocity of the wave.

Symbol:
    :code:`v`
"""

wave_angular_frequency = Symbol("wave_angular_frequency", angle_type / units.time)
r"""
Angular frequency of the wave.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

wave_amplitude = Symbol("wave_amplitude", units.length)
r"""
Amplitude of the wave.

Symbol:
    :code:`u_max`

Latex:
    :math:`u_\text{max}`
"""

law = Eq(
    wave_average_power,
    string_linear_density * wave_phase_velocity * wave_angular_frequency**2 * wave_amplitude**2 / 2,
)
r"""
:code:`P = 1/2 * mu * v * w^2 * u_max^2`

Latex:
    .. math::
        P = \frac{1}{2} \mu v \omega^2 u_\text{max}^2
"""


@validate_input(
    string_linear_density_=string_linear_density,
    wave_phase_velocity_=wave_phase_velocity,
    wave_angular_frequency_=wave_angular_frequency,
    wave_amplitude_=wave_amplitude,
)
@validate_output(wave_average_power)
def calculate_average_power(
    string_linear_density_: Quantity,
    wave_phase_velocity_: Quantity,
    wave_angular_frequency_: Quantity,
    wave_amplitude_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        string_linear_density: string_linear_density_,
        wave_phase_velocity: wave_phase_velocity_,
        wave_angular_frequency: wave_angular_frequency_,
        wave_amplitude: wave_amplitude_,
    })
    return Quantity(result)
