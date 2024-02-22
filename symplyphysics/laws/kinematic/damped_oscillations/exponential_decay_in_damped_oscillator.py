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
## For any positive value of the damping ratio, the position of a damped oscillator
## is subject to exponential decay, i.e. it is proportional to exp(-lambda*t) where
## lambda is the exponential decay constant.

# Law: lambda = omega * zeta
## lambda - exponential decay constant
## omega - undamped angular frequency of oscillator
## zeta - damping ratio

exponential_decay_constant = Symbol("exponential_decay_constant", 1 / units.time)
undamped_angular_frequency = Symbol("undamped_angular_frequency", angle_type / units.time)
damping_ratio = Symbol("damping_ratio", dimensionless)

law = Eq(exponential_decay_constant, undamped_angular_frequency * damping_ratio)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    undamped_angular_frequency_=undamped_angular_frequency,
    damping_ratio_=damping_ratio,
)
@validate_output(exponential_decay_constant)
def calculate_exponential_decay_constant(
    undamped_angular_frequency_: Quantity,
    damping_ratio_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        undamped_angular_frequency: undamped_angular_frequency_,
        damping_ratio: damping_ratio_,
    })
    return Quantity(result)
