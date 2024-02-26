from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The phase velocity of a wave is the rate at which the wave propagates in a medium.
## It is the velocity at which the phase of one frequency component of the wave travels.

# Law: v = lambda / T
## v - phase velocity of wave
## lambda - wavelength
## T - time period of wave

phase_velocity = Symbol("phase_velocity", units.length / units.time, positive=True)
wavelength = Symbol("wavelength", units.length, positive=True)
time_period = Symbol("time_period", units.time, positive=True)

law = Eq(phase_velocity, wavelength / time_period)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    wavelength_=wavelength,
    time_period_=time_period,
)
@validate_output(phase_velocity)
def calculate_phase_velocity(
    wavelength_: Quantity,
    time_period_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        wavelength: wavelength_,
        time_period: time_period_,
    })
    return Quantity(result)
