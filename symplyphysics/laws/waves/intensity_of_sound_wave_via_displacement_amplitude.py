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
## The intensity of a sound wave is the rate per unit area of energy transfer through or onto a surface.
## It is related to the displacement amplitude of the particles in the medium:

# Law: I = (1/2) * rho * v * w**2 * s_max**2
## I - intensity of sound wave
## rho - density of medium
## v - phase speed of sound wave
## w - angular frequency of sound wave
## s_max - displacement amplitude of particles in medium

intensity = Symbol("intensity", units.power / units.area)
density = Symbol("density", units.mass / units.volume)
phase_speed = Symbol("phase_speed", units.velocity)
angular_frequency = Symbol("angular_frequency", angle_type / units.time)
displacement_amplitude = Symbol("displacement_amplitude", units.length)

law = Eq(
    intensity,
    density * phase_speed * angular_frequency**2 * displacement_amplitude**2 / 2,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    density_=density,
    phase_speed_=phase_speed,
    angular_frequency_=angular_frequency,
    displacement_amplitude_=displacement_amplitude,
)
@validate_output(intensity)
def calculate_intensity(
    density_: Quantity,
    phase_speed_: Quantity,
    angular_frequency_: Quantity,
    displacement_amplitude_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        density: density_,
        phase_speed: phase_speed_,
        angular_frequency: angular_frequency_,
        displacement_amplitude: displacement_amplitude_,
    })
    return Quantity(result)
