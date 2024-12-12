from sympy import Eq, solve, sqrt
from sympy.physics.units import speed_of_light

from symplyphysics import (Quantity, Symbol, units, validate_input,
    validate_output)

# Description
# Time dilation is the difference in elapsed time as measured by two clocks,
# either due to a relative velocity between them (special relativity),
# or a difference in gravitational potential between their locations
# (general relativity).
# Law: t_rel = t / sqrt(1 - v**2 / c**2), where
# t_rel is time that has passed as measured by a stationary observer (relative time),
# t is time that has passed as measured by the traveling observer,
# v is velocity,
# c is speed of light.

# Links: formula in box <https://en.wikipedia.org/wiki/Time_dilation#Simple_inference>

moving_observer_time = Symbol("moving_observer_time", units.time)
velocity = Symbol("velocity", units.velocity)
relativistic_time = Symbol("relativistic_time", units.time)

law = Eq(relativistic_time, moving_observer_time / sqrt(1 - velocity**2 / speed_of_light**2))


@validate_input(moving_observer_time_=moving_observer_time, velocity_=velocity)
@validate_output(relativistic_time)
def calculate_relativistic_time(moving_observer_time_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, relativistic_time, dict=True)[0][relativistic_time]
    time_applied = result_expr.subs({
        moving_observer_time: moving_observer_time_,
        velocity: velocity_
    })
    return Quantity(time_applied)
