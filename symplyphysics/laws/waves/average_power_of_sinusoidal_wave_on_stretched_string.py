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
## The average power of a wave of any type is proportional to the square of its amplitude and to the
## square of its angular frequency.

# P_avg = mu * v * w**2 * u_max**2 / 2
## P_avg - average power, or rate of energy transmission, of wave
## mu - linear density of string
## v - phase velocity of wave
## w - angular frequency of wave
## u_max - amplitude of wave

wave_average_power = Symbol("wave_average_power", units.power)
string_linear_density = Symbol("string_linear_density", units.mass / units.length)
wave_phase_velocity = Symbol("wave_phase_velocity", units.velocity)
wave_angular_frequency = Symbol("wave_angular_frequency", angle_type / units.time)
wave_amplitude = Symbol("wave_amplitude", units.length)

law = Eq(
    wave_average_power,
    string_linear_density
    * wave_phase_velocity
    * wave_angular_frequency**2
    * wave_amplitude**2
    / 2,
)


def print_law() -> str:
    return print_expression(law)


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
