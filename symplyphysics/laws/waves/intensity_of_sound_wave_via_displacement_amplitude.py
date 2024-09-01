"""
Intensity of sound wave via displacement amplitude
==================================================

The intensity of a sound wave is the rate per unit area of energy transfer
through or onto a surface. It depends on the density of the medium, the phase
speed and the angular frequency of the wave and the amplitude of particles
in the medium.
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

wave_intensity = Symbol("wave_intensity", units.power / units.area)
"""
Intensity of the sound wave.

Symbol:
    :code:`I`
"""

medium_density = Symbol("medium_density", units.mass / units.volume)
r"""
Density of the medium in which the sound wave is being propagated.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

phase_speed = Symbol("phase_speed", units.velocity)
"""
Phase speed of the wave.

Symbol:
    :code:`v`
"""

angular_frequency = Symbol("angular_frequency", angle_type / units.time)
r"""
Angular frequency of the wave.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

displacement_amplitude = Symbol("displacement_amplitude", units.length)
r"""
Displacement amplitude of the particles in the medium.

Symbol:
    :code:`s_max`

Latex:
    :math:`s_\text{max}`
"""

law = Eq(
    wave_intensity,
    medium_density * phase_speed * angular_frequency**2 * displacement_amplitude**2 / 2,
)
r"""
:code:`I = (1 / 2) * rho * v * w^2 * s_max^2`

Latex:
    .. math::
        I = \frac{1}{2} \rho v \omega^2 s_\text{max}^2
"""


@validate_input(
    density_=medium_density,
    phase_speed_=phase_speed,
    angular_frequency_=angular_frequency,
    displacement_amplitude_=displacement_amplitude,
)
@validate_output(wave_intensity)
def calculate_intensity(
    density_: Quantity,
    phase_speed_: Quantity,
    angular_frequency_: Quantity,
    displacement_amplitude_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        medium_density: density_,
        phase_speed: phase_speed_,
        angular_frequency: angular_frequency_,
        displacement_amplitude: displacement_amplitude_,
    })
    return Quantity(result)
