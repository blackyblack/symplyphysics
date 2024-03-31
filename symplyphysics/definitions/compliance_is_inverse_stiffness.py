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
## Compliance, or flexibility, of a spring is the inverse of its stiffness and measures
## how flexible the spring is.

# Law: c = 1/k
## c - compliance of spring
## k - stiffness of spring

compliance = Symbol("compliance", units.length / units.force)
stiffness = Symbol("stiffness", units.force / units.length)

definition = Eq(compliance, 1 / stiffness)


def print_law() -> str:
    return print_expression(definition)


@validate_input(stiffness_=stiffness)
@validate_output(compliance)
def calculate_compliance(stiffness_: Quantity) -> Quantity:
    result = definition.rhs.subs(stiffness, stiffness_)
    return Quantity(result)
