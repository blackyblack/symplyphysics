from sympy import Eq, pi
from symplyphysics import (
    angle_type,
    dimensionless,
    Symbol,
    Quantity,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The interference of two waves is said to be fully constructive when the amplitude of the
## resulting wave is precisely the sum of the amplitudes of the comprising waves.

# Law: phi = 2 * pi * n
## phi - phase shift between interfering waves
## n - integer factor

phase_shift = Symbol("phase_shift", angle_type, real=True)
integer_factor = Symbol("integer_factor", dimensionless, integer=True)

law = Eq(phase_shift, 2 * pi * integer_factor)


def print_law() -> str:
    return print_expression(law)


@validate_input(integer_factor_=integer_factor)
@validate_output(phase_shift)
def calculate_phase_shift(integer_factor_: int) -> Quantity:
    result = law.rhs.subs(integer_factor, integer_factor_)
    return Quantity(result)
