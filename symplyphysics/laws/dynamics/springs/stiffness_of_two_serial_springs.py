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
## stiffness of the system of springs is equal to the sum of the stiffnesses of each spring.

# Law: k = k1 + k2
## k - total stiffness
## k1 - stiffness of first spring
## k2 - stiffness of second spring

# Condition
## - Springs must be Hookean, or linear-response, i.e. obey the Hooke's law.

total_stiffness = Symbol("total_stiffness", units.force / units.length)
first_stiffness = Symbol("first_stiffness", units.force / units.length)
second_stiffness = Symbol("second_stiffness", units.force / units.length)

law = Eq(total_stiffness, first_stiffness + second_stiffness)

# TODO: derive law from Hooke's law


def print_law() -> str:
    return print_expression(law)


@validate_input(
    first_stiffness_=first_stiffness,
    second_stiffness_=second_stiffness,
)
@validate_output(total_stiffness)
def calculate_total_stiffness(
    first_stiffness_: Quantity,
    second_stiffness_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        first_stiffness: first_stiffness_,
        second_stiffness: second_stiffness_,
    })
    return Quantity(result)
