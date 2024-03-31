from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    angle_type,
    dimensionless,
    validate_input,
    validate_output,
    print_expression,
)

# Description
## The damping ratio of an oscillator can be calculated via the exponential decay constant,
## which describes how fast the oscillations decay, and its undamped angular frequency.

# Law: zeta = lambda / omega
## zeta - damping ratio
## lambda - exponential decay constant
## omega - undamped angular frequency of oscillator

damping_ratio = Symbol("damping_ratio", dimensionless)
exponential_decay_constant = Symbol("exponential_decay_constant", 1 / units.time)
undamped_angular_frequency = Symbol("undamped_angular_frequency", angle_type / units.time)

law = Eq(damping_ratio, exponential_decay_constant / undamped_angular_frequency)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    exponential_decay_constant_=exponential_decay_constant,
    undamped_angular_frequency_=undamped_angular_frequency,
)
@validate_output(damping_ratio)
def calculate_damping_ratio(
    exponential_decay_constant_: Quantity,
    undamped_angular_frequency_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        exponential_decay_constant: exponential_decay_constant_,
        undamped_angular_frequency: undamped_angular_frequency_,
    })
    return Quantity(result)
