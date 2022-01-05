from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Newton's third law: Fr = -Fa
## Where:
## Fa - action force.
## Fr - reaction force.

force_action, force_reaction = symbols('force_action force_reaction')
law = Eq(force_reaction, -1 * force_action)

def print():
    return pretty(law, use_unicode=False)

@validate_input(force_action_=units.force)
@validate_output(units.force)
def calculate_force_reaction(force_action_: Quantity) -> Quantity:
    result_force_expr = solve(law, force_reaction, dict=True)[0][force_reaction]
    result_expr = result_force_expr.subs({force_action: force_action_})
    return expr_to_quantity(result_expr, 'force_reaction')
