from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Newton's third law: Fr = -Fa
## Where:
## Fa - action force.
## Fr - reaction force.

force_action = Symbol("force_action", units.force)
force_reaction = Symbol("force_reaction", units.force)

law = Eq(force_reaction, -1 * force_action)


def print_law() -> str:
    return print_expression(law)


@validate_input(force_action_=force_action)
@validate_output(force_reaction)
def calculate_force_reaction(force_action_: Quantity) -> Quantity:
    result_force_expr = solve(law, force_reaction, dict=True)[0][force_reaction]
    result_expr = result_force_expr.subs({force_action: force_action_})
    return Quantity(abs(result_expr))
