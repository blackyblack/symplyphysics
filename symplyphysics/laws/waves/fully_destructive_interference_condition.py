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
## The interference of two waves is said to be fully destructive when the amplitude of the
## resulting wave is zero, i.e. the two waves cancel each other out.

# Law: phi = (1 + 2 * n) * pi
## phi - phase shift between interfering waves
## n - integer factor

phase_shift = Symbol("phase_shift", angle_type, real=True)
integer_factor = Symbol("integer_factor", dimensionless, integer=True)

law = Eq(phase_shift, (1 + 2 * integer_factor) * pi)


def print_law() -> str:
    return print_expression(law)


@validate_input(integer_factor_=integer_factor)
@validate_output(phase_shift)
def calculate_phase_shift(integer_factor_: int) -> Quantity:
    result = law.rhs.subs(integer_factor, integer_factor_)
    return Quantity(result)
