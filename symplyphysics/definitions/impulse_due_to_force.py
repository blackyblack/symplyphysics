from sympy import Eq, Integral
from symplyphysics import (
    units,
    Quantity,
    Function,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The impulse due to force exerted on body during collision is the measure of both the
## magnitude and duration of the collision.

# Law: J = Integral(F(t), (t, t_i, t_f))
## J - magnitude of impulse due to force F
## F - magnitude of force
## t - time
## t_i, t_f - initial and final time of collision, respectively

impulse = Symbol("impulse", units.momentum)
force_function = Function("force_function", units.force)
time = Symbol("time", units.time)
time_start = Symbol("time_start", units.time)
time_end = Symbol("time_end", units.time)

law = Eq(impulse, Integral(force_function(time), (time, time_start, time_end)))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    force_start_=force_function,
    force_end_=force_function,
    time_start_=time_start,
    time_end_=time_end,
)
@validate_output(impulse)
def calculate_impulse(
    force_start_: Quantity,
    force_end_: Quantity,
    time_start_: Quantity,
    time_end_: Quantity,
) -> Quantity:
    force_function_ = force_start_ + (force_end_ - force_start_) / (time_end_ -
        time_start_) * (time - time_start_)
    result = law.rhs.subs({
        force_function(time): force_function_,
        time_start: time_start_,
        time_end: time_end_,
    }).doit()
    return Quantity(result)
