from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## If two springs are side-to-side to one another, i.e. connected in parallel, the total
## spring constant of the system of springs is equal to the sum of the spring constants
## of each spring.

# Law: k = k1 + k2
## k - total spring constant
## k1 - spring constant of first spring
## k2 - spring constant of second spring

# Condition
## - Springs must be Hookean, or linear-response, i.e. obey the Hooke's law.

total_spring_constant = Symbol("total_spring_constant", units.force / units.length)
first_spring_constant = Symbol("first_spring_constant", units.force / units.length)
second_spring_constant = Symbol("second_spring_constant", units.force / units.length)

law = Eq(total_spring_constant, first_spring_constant + second_spring_constant)

# TODO: derive law from Hooke's law


def print_law() -> str:
    return print_expression(law)


@validate_input(
    first_spring_constant_=first_spring_constant,
    second_spring_constant_=second_spring_constant,
)
@validate_output(total_spring_constant)
def calculate_spring_constant(
    first_spring_constant_: Quantity,
    second_spring_constant_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        first_spring_constant: first_spring_constant_,
        second_spring_constant: second_spring_constant_,
    })
    return Quantity(result)
