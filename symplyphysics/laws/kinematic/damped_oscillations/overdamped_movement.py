from sympy import Eq, exp
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    angle_type,
    validate_input,
    validate_output,
    print_expression,
)

# Description
## In the presence of a damping force in the oscillating system, the system's behaviour
## depends on the value of the damping ratio. When it is greater than 1, the system is unable
## to oscillate and instead slowly returns to equilibrium (exponentially decays to steady
## state). This behaviour is also known as overdamping.

# Law: x(t) = exp(-lambda*t) * (C1 * exp(|omega_damped|*t) + C2 * exp(-|omega_damped|*t))
## x(t) - position of overdamped oscillator
## t - time
## lambda - see [exponential decay constant](./exponential_decay_in_damped_oscillator.py)
## omega_damped - see [damped angular frequency](./damped_angular_frequency.py)

displacement = Function("displacement", units.length)
time = Symbol("time", units.time, positive=True)
exponential_decay_constant = Symbol("exponential_decay_constant", 1 / units.time)
damped_angular_frequency = Symbol("damped_angular_frequency", angle_type / units.time)
c1 = Symbol("C1", units.length)
c2 = Symbol("C2", units.length)

decay_factor = exp(-1 * exponential_decay_constant * time)
frequency_factor = exp(abs(damped_angular_frequency) * time)

law = Eq(
    displacement(time),
    decay_factor * (c1 * frequency_factor + c2 / frequency_factor)
)


def print_law() -> str:
    return print_expression(law)


# TODO: derive from damped oscillator equation directly


@validate_input(
    exponential_decay_constant_=exponential_decay_constant,
    damped_angular_frequency_=damped_angular_frequency,
    c1_=c1, c2_=c2, time_=time,
)
@validate_output(displacement)
def calculate_displacement(
    exponential_decay_constant_: Quantity,
    damped_angular_frequency_: Quantity,
    c1_: Quantity,
    c2_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        exponential_decay_constant: exponential_decay_constant_,
        damped_angular_frequency: damped_angular_frequency_,
        c1: c1_,
        c2: c2_,
        time: time_,
    })
    return Quantity(result)
