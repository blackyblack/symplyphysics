from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output, dimensionless)

# Description
## Since the hydraulic press is a mechanism, its operation can be characterized by a coefficient of efficiency.

## Law: n = (F2 * h2) / (F1 * h1)
## Where:
## F2 is the force spent on useful work
## h2 is the useful movement
## F1 is expended force
## h1 is expended height
## n is coefficient of efficiency

## Conditions
## The efficiency is calculated if we are talking about an imperfect hydraulic press.


useful_force = Symbol("useful_force", units.force)
useful_height = Symbol("useful_height", units.length)
expended_force = Symbol("expended_force", units.force)
expended_height = Symbol("expended_height", units.length)
efficiency = Symbol("efficiency", dimensionless)


law = Eq(efficiency, (useful_force * useful_height) / (expended_force * expended_height))


def print_law() -> str:
    return print_expression(law)


@validate_input(useful_force_=useful_force, useful_height_=useful_height, expended_force_=expended_force, expended_height_=expended_height)
@validate_output(efficiency)
def calculate_efficiency(useful_force_: Quantity, useful_height_, expended_force_, expended_height_: Quantity) -> Quantity:
    result_expr = solve(law, efficiency, dict=True)[0][efficiency]
    result_efficiency = result_expr.subs({
        useful_force: useful_force_,
        useful_height: useful_height_,
        expended_force: expended_force_,
        expended_height: expended_height_
    })

    return Quantity(result_efficiency)
