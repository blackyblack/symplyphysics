from sympy import Eq, exp, cos
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)

# Description
## In the presence of a damping force in the oscillating system, the system's behaviour
## depends on the value of the damping ratio. When it is less than 1, the system oscillates
## with a slightly different frequency than the undamped case with the amplitude decreasing
## to zero. This behaviour is also known as underdamping.

# Law: x(t) = x_max * exp(-lambda*t) * cos(omega_damped*t + phi)
## x(t) - position of underdamped oscillator
## t - time
## x_max - amplitude of oscillations
## lambda - see [exponential decay constant](./exponential_decay_in_damped_oscillator.py)
## omega_damped - see [damped angular frequency](./damped_angular_frequency.py)
## phi - phase lag of oscillations

displacement = Function("displacement", units.length, real=True)
time = Symbol("time", units.time, positive=True)
amplitude = Symbol("amplitude", units.length, positive=True)
exponential_decay_constant = Symbol("exponential_decay_constant", 1 / units.time, positive=True)
damped_angular_frequency = Symbol("damped_angular_frequency", angle_type / units.time, positive=True)
phase_lag = Symbol("phase_lag", angle_type, real=True)

decay_factor = exp(-1 * exponential_decay_constant * time)
frequency_factor = cos(damped_angular_frequency * time + phase_lag)

law = Eq(
    displacement(time),
    amplitude * decay_factor * frequency_factor,
)

# TODO: relate to [damped oscillator equation](../../../definitions/damped_harmonic_oscillator_equation.py)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    amplitude_=amplitude,
    exponential_decay_constant_=exponential_decay_constant,
    damped_angular_frequency_=damped_angular_frequency,
    phase_lag_=phase_lag,
    time_=time,
)
@validate_output(displacement)
def calculate_displacement(
    amplitude_: Quantity,
    exponential_decay_constant_: Quantity,
    damped_angular_frequency_: Quantity,
    phase_lag_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        amplitude: amplitude_,
        exponential_decay_constant: exponential_decay_constant_,
        damped_angular_frequency: damped_angular_frequency_,
        phase_lag: phase_lag_,
        time: time_,
    })
    return Quantity(result)
