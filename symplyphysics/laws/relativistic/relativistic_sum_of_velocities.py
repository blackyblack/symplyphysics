from sympy import Eq, solve
from sympy.physics.units import speed_of_light

from symplyphysics import (Quantity, Symbol, units, validate_input,
    validate_output)

# Description
# In relativistic mechanics, if a body moves relative to a moving reference
# system with velocity v1, and the velocity of the system relative to
# the observer is v2, then the velocity of the body relative
# to the observer will be:
# Law: v = (v1 + v2) / (1 + (v1 * v2) / c**2), where
# v1 is first velocity,
# v2 is second velocity,
# c is speed of light,
# v is relativistic sum of velocities.

# Links: Wikipedia <https://en.wikipedia.org/wiki/Velocity-addition_formula#Special_relativity>

first_velocity = Symbol("first_velocity", units.velocity)
second_velocity = Symbol("second_velocity", units.velocity)
resulting_velocity = Symbol("resulting_velocity", units.velocity)

law = Eq(resulting_velocity, (first_velocity + second_velocity) / (1 +
    (first_velocity * second_velocity) / speed_of_light**2))


@validate_input(
    first_velocity_=first_velocity,
    second_velocity_=second_velocity,
)
@validate_output(resulting_velocity)
def calculate_velocity(first_velocity_: Quantity, second_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, resulting_velocity)[0]
    velocity_applied = result_expr.subs({
        first_velocity: first_velocity_,
        second_velocity: second_velocity_,
    })
    return Quantity(velocity_applied)
