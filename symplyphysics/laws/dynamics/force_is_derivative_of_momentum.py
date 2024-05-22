from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)

# Description
## Newton's second law of motion can be generalized in terms of linear momentum.

# Law: dp/dt = F
## p - momentum projection
## F - force projection
## t - time

# Notes
## - Works in the relativistic case as well
## - Also see its [vector counterpart](./vector/force_is_derivative_of_momentum.py)

momentum = Function("momentum", units.momentum)
force = Function("force", units.force)
time = Symbol("time", units.time)

law = Eq(Derivative(momentum(time), time), force(time))


@validate_input(
    momentum_change_=momentum,
    time_change_=time,
)
@validate_output(force)
def calculate_force(
    momentum_change_: Quantity,
    time_change_: Quantity,
) -> Quantity:
    momentum_ = (momentum_change_ / time_change_) * time
    result = law.lhs.subs(momentum(time), momentum_).doit()
    return Quantity(result)
