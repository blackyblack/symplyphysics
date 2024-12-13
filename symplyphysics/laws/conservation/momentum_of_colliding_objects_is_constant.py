from sympy import (Derivative, Eq, dsolve)
from symplyphysics import (units, Quantity, Symbol, Function, print_expression, validate_input,
    validate_output)

# Description
## If there is no external force applied to system of objects, the summary momentum of this system remains constant
## during and after any interactions between objects.
## Summary momentum of the system is the sum of momentums of every object in this system.
## Also applicable for reactive engine simulation.

# Law: dP/dt = 0
## Where:
## P - summary momentum of colliding objects,
## t - time,
## dP/dt - momentum derivative over time.

# Conditions:
## 1. The system is closed

# Links: Wikipedia <https://en.wikipedia.org/wiki/Momentum#Conservation>

time = Symbol("time", units.time)
momentum = Function("momentum", units.momentum)

law = Eq(Derivative(momentum(time), time), 0)


@validate_input(momentum_before_=momentum)
@validate_output(momentum)
def calculate_momentum_after(momentum_before_: Quantity) -> Quantity:
    solved = dsolve(law, momentum(time))
    result_expr = solved.subs("C1", momentum_before_).rhs
    return Quantity(result_expr)
