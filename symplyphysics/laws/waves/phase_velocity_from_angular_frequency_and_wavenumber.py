from sympy import Eq, symbols, solve
from symplyphysics import (
    units,
    angle_type,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import angular_wavenumber_is_inverse_wavelength as wavenumber_def
from symplyphysics.laws.waves import wavelength_from_wave_speed_and_period as velocity_law
from symplyphysics.laws.kinematic import period_from_angular_frequency as frequency_law

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


# Derive from definition of phase velocity via wavelength and time period

_wavelength = solve(wavenumber_def.definition, wavenumber_def.wavelength)[0].subs(
    wavenumber_def.angular_wavenumber, angular_wavenumber
)

_time_period = solve(frequency_law.law, frequency_law.period)[0].subs(
    frequency_law.circular_frequency, angular_frequency
)

_phase_velocity = solve(velocity_law.law, velocity_law.propagation_speed)[0].subs({
    velocity_law.wavelength: _wavelength,
    velocity_law.oscillation_period: _time_period,
})

assert expr_equals(law.rhs, _phase_velocity)


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
