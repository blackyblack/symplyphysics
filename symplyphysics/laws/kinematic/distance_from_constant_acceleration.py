from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Motion with constant acceleration is described by the equation: S = S0 + V0 * t + (a * t^2)/2
## Where:
## S0 is the distance in the start moment
## V0 is the velocity in the start moment
## a is the acceleration.
## t is time.

initial_distance = Symbol("initial_distance", units.velocity)
time = Symbol("time", units.time)
acceleration = Symbol("acceleration", units.acceleration)
initial_velocity = Symbol("initial_velocity", units.velocity)
distance = Symbol("distance", units.length)

law = Eq(distance, initial_distance + initial_velocity * time + acceleration * time ** 2 / 2)


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_distance_=initial_distance, initial_velocity_=initial_velocity, acceleration_=acceleration, time_=time)
@validate_output(distance)
def calculate_distance(initial_distance_: Quantity, initial_velocity_: Quantity, acceleration_: Quantity, time_: Quantity) -> Quantity:
    solved = solve(law, distance, dict=True)[0][distance]
    result_expr = solved.subs({
        initial_distance: initial_distance_,
        initial_velocity: initial_velocity_,
        acceleration: acceleration_,
        time: time_
    })
    return Quantity(result_expr)
