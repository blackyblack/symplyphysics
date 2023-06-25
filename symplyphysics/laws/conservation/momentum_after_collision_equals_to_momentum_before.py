from sympy import (Eq, solve, dsolve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, Function, print_expression,
    validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.conservation import momentum_of_colliding_objects_is_constant as constant_momentum

# Description
## If there is no external force applied to system of objects, the summary momentum of this system remains constant
## during and after any interactions between objects.
## See [momentum_of_colliding_objects_is_constant](./momentum_of_colliding_objects_is_constant.py) for more information.

# Law: P(t1) = P(t0)
## Where:
## P - summary momentum of system of objects,
## t1 - point of time, when interaction of objects in a system occured,
## t0 - initial time.

time_before = Symbol("time_before", units.time)
time_after = Symbol("time_after", units.time)
momentum = Function("momentum", units.momentum)

law = Eq(momentum(time_after), momentum(time_before))

# Derive the same law from constant momentum

## dsolve() shows that solution is constant C1
dsolved = dsolve(constant_momentum.law, constant_momentum.momentum(constant_momentum.time))

energy_before_eq = dsolved.subs(constant_momentum.time, time_before)
energy_before_eq = energy_before_eq.subs(constant_momentum.momentum(time_before),
    momentum(time_before))
energy_after_eq = dsolved.subs(constant_momentum.time, time_after)
energy_after_eq = energy_after_eq.subs(constant_momentum.momentum(time_after), momentum(time_after))

## Show that when energy is constant, energy_before equals to energy_after
energy_after_solved = solve([energy_after_eq, energy_before_eq], (momentum(time_after), "C1"),
    dict=True)[0][momentum(time_after)]
assert expr_equals(energy_after_solved, law.rhs)


def print() -> str:
    return print_expression(law)


@validate_input(momentum_before_=momentum)
@validate_output(momentum)
def calculate_momentum_after(momentum_before_: Quantity) -> Quantity:
    solved = solve(law, momentum(time_after), dict=True)[0][momentum(time_after)]
    result_expr = solved.subs(momentum(time_before), momentum_before_)
    return expr_to_quantity(result_expr)
