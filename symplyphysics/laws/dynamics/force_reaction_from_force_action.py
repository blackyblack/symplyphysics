from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Newton's third law: Fr = -Fa
## Where:
## Fa - action force.
## Fr - reaction force.

force_action = Symbol("force_action", units.force)
force_reaction = Symbol("force_reaction", units.force)

law = Eq(force_reaction, -1 * force_action)

def print(expr: Expr) -> str:
    symbols = [force_action, force_reaction]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(force_action_=force_action)
@validate_output_symbol(force_reaction)
def calculate_force_reaction(force_action_: Quantity) -> Quantity:
    result_force_expr = solve(law, force_reaction, dict=True)[0][force_reaction]
    result_expr = result_force_expr.subs({force_action: force_action_})
    return expr_to_quantity(abs(result_expr))
