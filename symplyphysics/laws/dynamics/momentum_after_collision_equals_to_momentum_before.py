from sympy import (Eq, solve)
from symplyphysics import (
    units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol
)

# Description
## P_after = P_before
## Where:
## P_after - summary momentum of system of objects after interaction
## P_before - initial summary momentum.

## In other words, if there is no external force is applied to system of objects, the summary momentum of this system remains constant during and after any interactions between objects
## Summary momentum of the system is the sum of momentums of every member of this system.
## Also applicable for reactive engine simulation

momentum_before = Symbol("momentum_before", units.momentum)
momentum_after = Symbol("momentum_after", units.momentum)

law = Eq(momentum_after, momentum_before)

def print() -> str:
    return print_expression(law)

@validate_input_symbols(momentum_before_=momentum_before)
@validate_output_symbol(momentum_after)
def calculate_momentum_after(momentum_before_: Quantity) -> Quantity:
    solved = solve(law, momentum_after, dict=True)[0][momentum_after]    
    result_expr = solved.subs(momentum_before, momentum_before_)
    return expr_to_quantity(result_expr)
