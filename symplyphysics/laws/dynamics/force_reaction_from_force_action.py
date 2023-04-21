from sympy import (Eq, solve)
from symplyphysics import (
    units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol
)

# Description
## Newton's third law: Fr = -Fa
## Where:
## Fa - action force.
## Fr - reaction force.

force_action = Symbol("force_action", units.force)
force_reaction = Symbol("force_reaction", units.force)

law = Eq(force_reaction, -1 * force_action)

def print() -> str:
    return print_expression(law)

@validate_input_symbols(force_action_=force_action)
@validate_output_symbol(force_reaction)
def calculate_force_reaction(force_action_: Quantity) -> Quantity:
    result_force_expr = solve(law, force_reaction, dict=True)[0][force_reaction]
    result_expr = result_force_expr.subs({force_action: force_action_})
    return expr_to_quantity(abs(result_expr))
