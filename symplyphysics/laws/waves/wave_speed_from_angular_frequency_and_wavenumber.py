from sympy import Eq
from symplyphysics import (
    units,
    angle_type,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The speed of a traveling wave can be found using its angular frequency and wavenumber.

# Law: v = w/k
## v - wave speed
## w - angular frequency of wave
## k - angular wavenumber of wave

wave_speed = Symbol("wave_speed", units.length / units.time, real=True)
angular_frequency = Symbol("angular_frequency", angle_type / units.time, positive=True)
angular_wavenumber = Symbol("angular_wavenumber", angle_type / units.length, positive=True)

law = Eq(wave_speed, angular_frequency / angular_wavenumber)


# TODO: derive from wave phase equation


def print_law() -> str:
    return print_expression(law)


@validate_input(
    angular_frequency_=angular_frequency,
    angular_wavenumber_=angular_wavenumber,
)
@validate_output(wave_speed)
def calculate_wave_speed(
    angular_frequency_: Quantity,
    angular_wavenumber_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        angular_frequency: angular_frequency_,
        angular_wavenumber: angular_wavenumber_,
    })
    return Quantity(result)
