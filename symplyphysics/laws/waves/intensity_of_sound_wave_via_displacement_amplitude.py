"""
Intensity of sound wave via displacement amplitude
==================================================

The intensity of a sound wave is the rate per unit area of energy transfer
through or onto a surface. It depends on the density of the medium, the phase
speed and the angular frequency of the wave and the amplitude of particles
in the medium.

**Links:**

#. Equation 17-27 on p. 489 of "Fundamentals of Physics" by David Halladay et al., 10th Ed.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

wave_intensity = symbols.intensity
"""
:symbols:`intensity` of the sound wave.
"""

medium_density = symbols.density
"""
:symbols:`density` of the medium in which the sound wave is being propagated.
"""

phase_speed = symbols.phase_speed
"""
:symbols:`phase_speed` of the wave.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the wave.
"""

displacement_amplitude = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="s_max",
    display_latex="s_\\text{max}",
)
"""
Displacement amplitude of the particles in the medium. See :symbols:`euclidean_distance`.
"""

law = Eq(
    wave_intensity,
    medium_density * phase_speed * angular_frequency**2 * displacement_amplitude**2 / 2,
)
"""
:laws:symbol::

:laws:latex::
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
