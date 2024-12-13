"""
Pressure amplitude in sound wave
================================

Sound waves cause a pressure change of the medium from the equilibrium pressure.
This change is proportional to the speed of sound and the density of the medium,
the angular frequency of the wave and the displacement of particles in the medium.

**Links:**

#. Equation 17-14 on p. 484 of "Fundamentals of Physics" by David Halladay et al., 10th Ed.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

pressure_amplitude = clone_as_symbol(
    symbols.pressure,
    display_symbol="Delta(p)_max",
    display_latex="(\\Delta p)_\\text{max}",
)
"""
Amplitude of :symbols:`pressure` change.
"""

speed_of_sound = symbols.speed
"""
:symbols:`speed` of sound in the medium.
"""

medium_density = symbols.density
"""
:symbols:`density` of the medium.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the sound wave.
"""

displacement_amplitude = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="s_max",
    display_latex="s_\\text{max}",
)
"""
Displacement amplitude of particles in the medium. See :symbols:`euclidean_distance`.
"""

law = Eq(
    pressure_amplitude,
    speed_of_sound * medium_density * angular_frequency * displacement_amplitude,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    speed_of_sound_=speed_of_sound,
    density_of_medium_=medium_density,
    angular_frequency_=angular_frequency,
    displacement_amplitude_=displacement_amplitude,
)
@validate_output(pressure_amplitude)
def calculate_pressure_amplitude(
    speed_of_sound_: Quantity,
    density_of_medium_: Quantity,
    angular_frequency_: Quantity,
    displacement_amplitude_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        speed_of_sound: speed_of_sound_,
        medium_density: density_of_medium_,
        angular_frequency: angular_frequency_,
        displacement_amplitude: displacement_amplitude_,
    })
    return Quantity(result)
