from sympy import Eq, sin, cos, symbols, Function as SymFunction
from symplyphysics import (
    units,
    angle_type,
    Symbol,
    Quantity,
    print_expression,
    validate_input,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.quantity_decorator import validate_output_same
from symplyphysics.laws.waves import phase_of_traveling_wave as phase_law

# Description
## A standing, or stationary, wave is the result of the interference of two identical waves
## moving in the opposite direction. For a wave occuring along a string, the standing wave
## is described by the following equation:

# Law: u(x, t) = 2 * u_max * sin(k * x) * cos(w * t)
## u(x, t) - displacement from rest (position, pressure, electric field, etc)
## u_max - wave amplitude
## k - angular wavenumber
## x - spatial variable
## w - angular frequency
## t - time

# Note
## - In this law we assume that the standing wave is composed of two identical sinusoidal traveling
##   waves whose form is described by the expression `u_max * sin(k*x +/- w*t)`
## - This is no longer a traveling wave because for that it should be moving in one direction

displacement = symbols("displacement", cls=SymFunction, real=True)
amplitude = symbols("amplitude", positive=True)
angular_wavenumber = Symbol("angular_wavenumber", angle_type / units.length, positive=True)
position = Symbol("position", units.length, real=True)
angular_frequency = Symbol("angular_frequency", angle_type / units.time, positive=True)
time = Symbol("time", units.time, real=True)

law = Eq(displacement(position, time),
    (2 * amplitude) * sin(angular_wavenumber * position) * cos(angular_frequency * time))

# Derive from the sum of two traveling waves

_wave_phase = phase_law.law.rhs.subs({
    phase_law.angular_wavenumber: angular_wavenumber,
    phase_law.position: position,
    phase_law.angular_frequency: angular_frequency,
    phase_law.time: time,
})

# The form of the waves is taken from the note above
_forward_wave = amplitude * sin(_wave_phase)

_backward_wave = _forward_wave.subs(angular_frequency, -1 * angular_frequency)

_sum_of_waves = (_forward_wave + _backward_wave).simplify()

assert expr_equals(_sum_of_waves, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    angular_wavenumber_=angular_wavenumber,
    position_=position,
    angular_frequency_=angular_frequency,
    time_=time,
)
@validate_output_same("amplitude_")
def calculate_displacement(
    amplitude_: Quantity,
    angular_wavenumber_: Quantity,
    position_: Quantity,
    angular_frequency_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        amplitude: amplitude_,
        angular_wavenumber: angular_wavenumber_,
        position: position_,
        angular_frequency: angular_frequency_,
        time: time_,
    })
    return Quantity(result)
