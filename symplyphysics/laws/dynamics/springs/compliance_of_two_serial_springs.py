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
## If two springs are connected end-to-end to one another, i.e. connected in series,
## the total compliance of the system of springs is the sum of the compliances of
## each spring.

# Law: c = c1 + c2
## c - total compliance
## c1 - compliance of first spring
## c2 - compliance of second spring

# Condition
## - Springs must be Hookean, or linear-response, i.e. obey the Hooke's law.

total_compliance = Symbol("total_compliance", units.length / units.force)
first_compliance = Symbol("first_compliance", units.length / units.force)
second_compliance = Symbol("second_compliance", units.length / units.force)

law = Eq(total_compliance, first_compliance + second_compliance)

# TODO: derive law from Hooke's law


def print_law() -> str:
    return print_expression(law)


@validate_input(
    first_compliance_=first_compliance,
    second_compliance_=second_compliance,
)
def calculate_compliance(
    first_compliance_: Quantity,
    second_compliance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        first_compliance: first_compliance_,
        second_compliance: second_compliance_,
    })
    return Quantity(result)
