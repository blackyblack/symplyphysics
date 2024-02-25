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
## The phase velocity of a wave is the rate at which the wave propagates in a medium.
## It is the velocity at which the phase of one frequency component of the wave travels.

# Law: v = w / k
## v - phase velocity of wave
## w - angular frequency of wave
## k - angular wavenumber of wave

phase_velocity = Symbol("phase_velocity", units.length / units.time, real=True)
angular_frequency = Symbol("angular_frequency", angle_type / units.time, positive=True)
angular_wavenumber = Symbol("angular_wavenumber", angle_type / units.length, positive=True)

law = Eq(phase_velocity, angular_frequency / angular_wavenumber)


# TODO: derive from wave phase equation


def print_law() -> str:
    return print_expression(law)


@validate_input(
    angular_frequency_=angular_frequency,
    angular_wavenumber_=angular_wavenumber,
)
@validate_output(phase_velocity)
def calculate_phase_velocity(
    angular_frequency_: Quantity,
    angular_wavenumber_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        angular_frequency: angular_frequency_,
        angular_wavenumber: angular_wavenumber_,
    })
    return Quantity(result)
