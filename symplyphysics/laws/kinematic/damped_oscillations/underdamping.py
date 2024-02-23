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
from symplyphysics.core.symbols.quantities import scale_factor

# Description
## In the presence of a damping force in the oscillating system, the system's behaviour
## depends on the value of the damping ratio. When it is less than 1, the system oscillates
## with a slightly different frequency than in the undamped case, and its amplitude decreasing
## to zero. This behaviour is also known as underdamping.

# Law: x(t) = x_max * exp(-lambda*t) * cos(omega_damped*t + phi)
## x(t) - position of damped oscillator
## x_max - amplitude of damped oscillations
## lambda - exponential decay constant
## omega_damped - damped angular frequency
## phi - phase lag

# Conditions
## - System is underdamped, i.e. its damping ratio is less than 1.

displacement = Function("displacement", units.length, real=True)
time = Symbol("time", units.time, nonnegative=True)
amplitude = Symbol("amplitude", units.length, positive=True)
exponential_decay_constant = Symbol("exponential_decay_constant", 1 / units.time, positive=True)
undamped_angular_frequency = Symbol("undamped_angular_frequency", angle_type / units.time, positive=True)
phase_lag = Symbol("phase_lag", angle_type, real=True)

law = Eq(
    displacement(time),
    amplitude
    * exp(-1 * exponential_decay_constant * time)
    * cos(undamped_angular_frequency * time + phase_lag)
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    amplitude_=amplitude,
    exponential_decay_constant_=exponential_decay_constant,
    undamped_angular_frequency_=undamped_angular_frequency,
    phase_lag_=phase_lag,
    time_=time,
)
@validate_output(displacement)
def calculate_displacement(
    amplitude_: Quantity,
    exponential_decay_constant_: Quantity,
    undamped_angular_frequency_: Quantity,
    phase_lag_: Quantity | float,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        amplitude: amplitude_,
        exponential_decay_constant: exponential_decay_constant_,
        undamped_angular_frequency: undamped_angular_frequency_,
        phase_lag: scale_factor(phase_lag_),
        time: time_,
    })
    return Quantity(result)
