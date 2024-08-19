from sympy import Eq, sqrt
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
## A damped oscillator performs oscillations with a frequency proportional
## to that of an undamped oscillator, i.e. to the frequency of an oscillator when
## there is no damping force acting on the system.

# Law: omega_damped = omega * sqrt(1 - zeta**2)
## omega_damped - angular frequency of damped oscillator
## omega - undamped angular frequency
## zeta - damping ratio

# Note
## In the case of overdamping when the damping ratio is greater than 1, the value of damped
## angular frequency becomes imaginary, which indicates the absense of oscillations.

damped_angular_frequency = Symbol("damped_angular_frequency", angle_type / units.time)
undamped_angular_frequency = Symbol("undamped_angular_frequency", angle_type / units.time)
damping_ratio = Symbol("damping_ratio", dimensionless)

law = Eq(damped_angular_frequency, undamped_angular_frequency * sqrt(1 - damping_ratio**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    undamped_angular_frequency_=undamped_angular_frequency,
    damping_ratio_=damping_ratio,
)
@validate_output(damped_angular_frequency)
def calculate_damped_angular_frequency(
    undamped_angular_frequency_: Quantity,
    damping_ratio_: float,
) -> Quantity:
    result = law.rhs.subs({
        undamped_angular_frequency: undamped_angular_frequency_,
        damping_ratio: damping_ratio_,
    })
    return Quantity(result)
