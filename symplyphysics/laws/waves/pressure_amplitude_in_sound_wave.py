from sympy import Eq
from symplyphysics import (
    units,
    angle_type,
    Symbol,
    Quantity,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Sound waves cause a pressure change of the medium from the equilibrium pressure.

# Law: delta_p_max = v * rho * w * s_max
## delta_p_max - amplitude of pressure change
## v - speed of sound in medium
## w - angular frequency of the sound wave
## s_max - displacement amplitude of medium particles

pressure_amplitude = Symbol("pressure_amplitude", units.pressure)
speed_of_sound = Symbol("speed_of_sound", units.velocity)
density_of_medium = Symbol("density_of_medium", units.mass / units.volume)
angular_frequency = Symbol("angular_frequency", angle_type / units.time)
displacement_amplitude = Symbol("displacement_amplitude", units.length)

law = Eq(
    pressure_amplitude,
    speed_of_sound * density_of_medium * angular_frequency * displacement_amplitude,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    speed_of_sound_=speed_of_sound,
    density_of_medium_=density_of_medium,
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
        density_of_medium: density_of_medium_,
        angular_frequency: angular_frequency_,
        displacement_amplitude: displacement_amplitude_,
    })
    return Quantity(result)
