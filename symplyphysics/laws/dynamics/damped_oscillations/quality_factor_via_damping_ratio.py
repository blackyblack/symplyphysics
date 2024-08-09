from sympy import Eq
from symplyphysics import (
    convert_to_float,
    dimensionless,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## In case of a damped oscillating system where the damping force is linearly proportional
## to the oscillator's frequency, the Q factor of the oscillator is inversely proportional
## to the system's damping ratio. Hence greater values of the Q factor correspond to lower
## values of the damping ratio and to lower damping.

# Law: Q = 1/(2 * zeta)
## Q - Q factor
## zeta - [damping ratio](../../../definitions/damped_harmonic_oscillator_equation.py)

quality_factor = Symbol("quality_factor", dimensionless)
damping_ratio = Symbol("damping_ratio", dimensionless)

law = Eq(quality_factor, 1 / (2 * damping_ratio))


def print_law() -> str:
    return print_expression(law)


@validate_input(damping_ratio_=damping_ratio)
@validate_output(quality_factor)
def calculate_quality_factor(damping_ratio_: float) -> float:
    result = law.rhs.subs(damping_ratio, damping_ratio_)
    return convert_to_float(result)
