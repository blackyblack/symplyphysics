from symplyphysics import (
    symbols, Function, Eq, pretty, dsolve, solve, units, Quantity, validate_input, validate_output, expr_to_quantity, simplify
)

from symplyphysics.definitions import velocity_is_movement_derivative as velocity_definition
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration_definition

# Description
## Accelerated movement is the kind of movement when object has constant acceleration (e.g with the constant force applied to object).
## The movement is S = V0*t + a * t**2/2, where
## S is distance in the moment of time t,
## V0 is the velocity in the moment of time 0 (initial velocity),
## a is acceleration,
## and t is time.

movement_time = symbols('movement_time')
constant_acceleration = symbols('constant_acceleration', constant=True)
initial_velocity = symbols('initial_velocity')
distance_function = symbols('distance_function', cls=Function)

law = Eq(distance_function(movement_time), initial_velocity * movement_time + constant_acceleration * movement_time**2 / 2)

# Derive the same law from velocity and acceleration definitions

constant_acceleration_definition = acceleration_definition.definition.subs({acceleration_definition.acceleration(acceleration_definition.time): constant_acceleration, acceleration_definition.time: movement_time})
dsolved_velocity = dsolve(constant_acceleration_definition, acceleration_definition.velocity_function(movement_time))
constant_accelerated_velocity_function = dsolved_velocity.rhs

constant_accelerated_movement_definition = velocity_definition.definition.subs({velocity_definition.velocity_function(velocity_definition.moving_time): constant_accelerated_velocity_function, velocity_definition.moving_time: movement_time})
dsolved_movement = dsolve(constant_accelerated_movement_definition, velocity_definition.movement_function(movement_time))
constant_accelerated_movement_function = dsolved_movement.rhs

derived_law = Eq(distance_function(movement_time), constant_accelerated_movement_function)

# Prove that constant_accelerated_movement_function equals to law.rhs, given C1 = initial_velocity, C2 = 0
difference = simplify(derived_law.rhs.subs({'C1': initial_velocity, 'C2': 0}) - law.rhs)
assert(difference == 0)

def print():
    return pretty(law, use_unicode=False)

@validate_input(initial_velocity_=units.velocity, acceleration_=units.acceleration, time_=units.time)
@validate_output(units.length)
def calculate_distance(initial_velocity_: Quantity, acceleration_: Quantity, time_: Quantity) -> Quantity:        
    result_expr = solve(law, distance_function(movement_time), dict=True)[0][distance_function(movement_time)]
    result_expr_substituted = result_expr.subs({initial_velocity: initial_velocity_, constant_acceleration: acceleration_, movement_time: time_})
    return expr_to_quantity(result_expr_substituted, 'distance')
