from sympy import Eq, log
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The sound level of a sound wave is another way to describe the wave's intensity:

# Definition: beta = beta0 * log10(I / I0)
## beta - sound level
## beta0 = 10 dB - dimension factor in decibel, the unit of sound level
## I - intensity of sound wave
## I0 = 1e-12 W/m**2 - reference intensity level

sound_level = Symbol("sound_level", dimensionless)
intensity = Symbol("intensity", units.power / units.area)

reference_sound_level = Quantity(10)
reference_intensity = Quantity(1e-12 * units.watt / units.meter**2)

definition = Eq(
    sound_level,
    reference_sound_level * log(intensity / reference_intensity, 10)
)


def print_law() -> str:
    return print_expression(definition)


@validate_input(intensity_=intensity)
@validate_output(sound_level)
def calculate_sound_level(intensity_: Quantity) -> float:
    result = definition.rhs.subs(intensity, intensity_)
    return Quantity(result).scale_factor
