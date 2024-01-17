from sympy import Eq, solve, Derivative
from symplyphysics import (
    Symbol,
    units,
    dimensionless,
    print_expression,
    Quantity,
    validate_input,
    validate_output,
)

# Description
## Electric charge is quantized, i.e. restricted to certain values.

# Law: q = n*e
## q - electric charge
## n - an integer factor (from the set of whole numbers, i.e. positive, zero, or negative)
## e - elementary charge

charge = Symbol("charge", units.charge)
integer_factor = Symbol("factor", dimensionless)  # positive, zero or negative

law = Eq(charge, integer_factor * units.elementary_charge)


def print_law() -> str:
    return print_expression(law)


@validate_input(integer_factor_=integer_factor)
@validate_output(charge)
def calculate_charge(integer_factor_: int) -> Quantity:
    charge_expr = law.rhs
    charge_value = charge_expr.subs(integer_factor, integer_factor_)
    return Quantity(charge_value)
