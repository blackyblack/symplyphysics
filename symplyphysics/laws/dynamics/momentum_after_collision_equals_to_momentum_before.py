from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## P_after = P_before
## Where:
## P_after - summary momentum of system of objects after interaction
## P_before - initial summary momentum.

## In other words, if there is no external force is applied to system of objects, the summary momentum of this system remains constant during and after any interactions between objects
## Summary momentum of the system is the sum of momentums of every member of this system.
## Also applicable for reactive engine simulation

momentum_before, momentum_after = symbols('momentum_before momentum_after')
law = Eq(momentum_after, momentum_before)

def print():
    return pretty(law, use_unicode=False)

@validate_input(momentum_before_ = units.momentum)
@validate_output(units.momentum)
def calculate_momentum_after(momentum_before_: Quantity) -> Quantity:
    solved = solve(law, momentum_after, dict=True)[0][momentum_after]    
    result_expr = solved.subs({momentum_before : momentum_before_})
    return expr_to_quantity(result_expr, 'momentum_after')
