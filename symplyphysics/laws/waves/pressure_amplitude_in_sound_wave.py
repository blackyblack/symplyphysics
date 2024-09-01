"""
Pressure amplitude in sound wave
================================

Sound waves cause a pressure change of the medium from the equilibrium pressure.
This change is proportional to the speed of sound and the density of the medium,
the angular frequency of the wave and the displacement of particles in the medium.
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

pressure_amplitude = Symbol("pressure_amplitude", units.pressure)
r"""
Amplitude of pressure change.

Symbol:
    :code:`Delta(p)_max`

Latex:
    :math:`(\Delta p)_\text{max}`
"""

speed_of_sound = Symbol("speed_of_sound", units.velocity)
"""
Speed of sound in the medium.

Symbol:
    :code:`v`
"""

medium_density = Symbol("medium_density", units.mass / units.volume)
r"""
Density of the medium.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

angular_frequency = Symbol("angular_frequency", angle_type / units.time)
r"""
Angular frequency of the sound wave.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

displacement_amplitude = Symbol("displacement_amplitude", units.length)
r"""
Displacement amplitude of particles in the medium.

Symbol:
    :code:`s_max`

Latex:
    :math:`s_\text{max}`
"""

law = Eq(
    pressure_amplitude,
    speed_of_sound * medium_density * angular_frequency * displacement_amplitude,
)
r"""
:code:`Delta(p)_max = v * rho * w * s_max`

Latex:
    .. math::
        (\Delta p)_\text{max} = v \rho \omega s_\text{max}
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
