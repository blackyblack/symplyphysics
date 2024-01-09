from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Work is measured result of force applied. Mechanical work is the only reason for the object energy to be changed.
## Work is scalar value equal to force multiplied by movement.
## Use vector form of this law for non-collinear vectors of force and movement.
## Law: A = F * S, where
## A is mechanical work
## F is magnitude of force applied to object
## S is magnitude of movement caused by this force

# Conditions:
## Force and work vectors are collinear

work = Symbol("work", units.energy)
force = Symbol("force", units.force)
distance = Symbol("distance", units.length)

law = Eq(work, force * distance)


def print_law() -> str:
    return print_expression(law)


@validate_input(force_=force, distance_=distance)
@validate_output(work)
def calculate_work(force_: Quantity, distance_: Quantity) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    result_expr = result_work_expr.subs({force: force_, distance: distance_})
    return Quantity(result_expr)
