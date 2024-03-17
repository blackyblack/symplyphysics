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
## Any function `h(x, t)` that depends on position `x` and time `t` can describe a traveling
## wave if `x` and `t` appear solely in the form of `k*x + omega*t`, i.e. the function depends
## on the wave phase described below.

# Law: phi(x, t) = k * x - omega * t
## phi - wave phase
## k - angular wavenumber
## x - position
## omega - angular frequency
## t - time

# Note
## `omega` can be both positive (the wave travels in the positive direction of x-axis)
## or negative (the wave travels in the negative direction of x-axis)

# Condition
## - This is the case of a wave traveling in a single spatial dimension.
## - A constant phase lag is not taken into account.

wave_phase = Symbol("wave_phase", angle_type, real=True)
angular_wavenumber = Symbol("angular_wavenumber", angle_type / units.length, positive=True)
position = Symbol("position", units.length, real=True)
angular_frequency = Symbol("angular_frequency", angle_type / units.time, real=True)
time = Symbol("time", units.time, positive=True)

law = Eq(wave_phase, angular_wavenumber * position - angular_frequency * time)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    angular_wavenumber_=angular_wavenumber,
    position_=position,
    angular_frequency_=angular_frequency,
    time_=time,
)
@validate_output(wave_phase)
def calculate_wave_phase(
    angular_wavenumber_: Quantity,
    position_: Quantity,
    angular_frequency_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        angular_wavenumber: angular_wavenumber_,
        position: position_,
        angular_frequency: angular_frequency_,
        time: time_,
    })
    return Quantity(result)
