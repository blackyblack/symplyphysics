from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, Function, print_expression,
    validate_input_symbols, validate_output_symbol)

# Description
## If there is no external force is applied to system of objects, the summary momentum of this system remains constant
## during and after any interactions between objects.
## Summary momentum of the system is the sum of momentums of every object in this system.
## Also applicable for reactive engine simulation

# Law: P(t1) = P(t0)
## Where:
## P - summary momentum of system of objects,
## t1 - point of time, when interaction of objects in a system occured,
## t0 - initial time.

time_before = Symbol("time_before", units.time)
time_after = Symbol("time_after", units.time)
momentum = Function("momentum", units.momentum)

law = Eq(momentum(time_after), momentum(time_before))


def print() -> str:
    return print_expression(law)


@validate_input_symbols(momentum_before_=momentum)
@validate_output_symbol(momentum)
def calculate_momentum_after(momentum_before_: Quantity) -> Quantity:
    solved = solve(law, momentum(time_after), dict=True)[0][momentum(time_after)]
    result_expr = solved.subs(momentum(time_before), momentum_before_)
    return expr_to_quantity(result_expr)
