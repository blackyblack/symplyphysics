from sympy import Eq, sin, cos, symbols, Function as SymFunction
from symplyphysics import (
    units,
    angle_type,
    Symbol,
    Quantity,
    print_expression,
    validate_input,
)
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.core.quantity_decorator import validate_output_same

# Description
## If two waves are traveling in the same direction and have the same amplitude, period, and
## wavelength (and hence the same frequency and wavenumber), but differ in phase constant, the result
## is a single wave with the same period and wavelength, but its amplitude depends on the phase
## shift between the waves. If the shift is zero (or a multiple of 2*pi), the waves are exactly
## in phase and their interference is fully constructive; if it is pi (or pi with a multiple of 2*pi),
## they are exactly out of phase and their interference is fully destructive.

# Law: u_sum(x, t) = 2 * u_max * cos(phi / 2) * sin(k * x - w * t + phi / 2)
## u_sum - displacement of the resulting wave
## u_max - amplitude of the interfering waves
## phi - phase shift between the interfering waves
## k - [angular wavenumber](../../definitions/angular_wavenumber_is_inverse_wavelength.py)
## w - angular frequency
## x - position
## t - time

# Conditions
## - The waves are traveling in the same (or similar) directions
## - They have the same amplitude, angular wavenumber and angular frequency

# Note
## - The form of the first wave is `u_max * sin(k * x - w *t)` and of the second wave is
##   `u_max * sin(k * x - w * t + phi)`, i.e. the second wave is shifted by `phi` relative
##   to the first one.
## - The travel of the waves is unaffected by their interference.

displacement = symbols("displacement", cls=SymFunction, real=True)
amplitude = symbols("amplitude", positive=True)
phase_shift = Symbol("phase_shift", angle_type, real=True)
angular_wavenumber = Symbol("angular_wavenumber", angle_type / units.length, positive=True)
angular_frequency = Symbol("angular_frequency", angle_type / units.time, positive=True)
position = Symbol("position", units.length, real=True)
time = Symbol("time", units.time, real=True)

law = Eq(displacement(position, time), (2 * amplitude) * cos(phase_shift / 2) *
    sin(angular_wavenumber * position - angular_frequency * time + phase_shift / 2))

# TODO: derive from the sum of two waves


def print_law() -> str:
    return print_expression(law)


@validate_input(
    phase_shift_=phase_shift,
    angular_wavenumber_=angular_wavenumber,
    angular_frequency_=angular_frequency,
    position_=position,
    time_=time,
)
@validate_output_same("amplitude_")
def calculate_displacement(
    amplitude_: Quantity,
    phase_shift_: Quantity | float,
    angular_wavenumber_: Quantity,
    angular_frequency_: Quantity,
    position_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        amplitude: amplitude_,
        phase_shift: scale_factor(phase_shift_),
        angular_wavenumber: angular_wavenumber_,
        angular_frequency: angular_frequency_,
        position: position_,
        time: time_,
    })
    return Quantity(result)
