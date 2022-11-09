from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Accelerated velocity is time dependent and increases with time if acceleration is co-directed with velocity and decreases if they are counter-directed.
## Initial velocity is velocity at the start of observation, when time is 0.
## In other words, velocity V = V0 + at, where
## V is the velocity in the moment of time t
## V0 is the initial velocity (in the moment of time 0)
## a is the acceleration.

velocity = symbols('velocity')
time = symbols('time')
acceleration = symbols('acceleration')
initial_velocity = symbols('initial_velocity')

law = Eq(velocity, initial_velocity + acceleration * time)

def print():
    return pretty(law, use_unicode=False)

@validate_input(initial_velocity_ = units.velocity, acceleration_ = units.acceleration, time_ = units.time)
@validate_output(units.velocity)
def calculate_velocity(initial_velocity_: Quantity, acceleration_: Quantity, time_: Quantity) -> Quantity:
    result_velocity_expression = solve(law, velocity, dict=True)[0][velocity]
    result_expr = result_velocity_expression.subs({initial_velocity: initial_velocity_, acceleration: acceleration_, time: time_})
    return expr_to_quantity(result_expr, 'accelerated_velocity')
