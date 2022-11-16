from symplyphysics import (
    symbols, Function, Eq, pretty, dsolve
)

from symplyphysics.definitions import velocity_is_movement_derivative as velocity_definition
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration_definition

# Description
## Accelerated movement is the kind of movement when object has constant acceleration (e.g with the constant force applied to object).
## For example, free fall is accelerated movement. The acceleration in each moment of time is g.
## According to velocity and acceleration definitions, velocity is movement derivative, acceleration is velocity derivative.
## The law is: S = V0*t + a * t**2/2, where
## S is distance in the moment of time t
## V0 is the velocity in the moment of time 0 (initial velocity)
## a is acceleration
## and t is time.

moving_time = symbols('moving_time')
constant_acceleration = symbols('constant_acceleration', constant=True)
initial_velocity = symbols('initial_velocity')
distance_function = symbols('distance_function', cls=Function)

law = Eq(distance_function(moving_time), initial_velocity * moving_time + constant_acceleration * moving_time**2 / 2)

# Derive the same law from velocity and acceleration definitions

constant_acceleration_definition = acceleration_definition.definition.subs({acceleration_definition.acceleration(acceleration_definition.time): constant_acceleration, acceleration_definition.time: moving_time})
dsolved_velocity = dsolve(constant_acceleration_definition, acceleration_definition.velocity_function(moving_time))
constant_accelerated_velocity_function = dsolved_velocity.rhs

constant_accelerated_movement_definition = velocity_definition.definition.subs({velocity_definition.velocity_function(velocity_definition.moving_time): constant_accelerated_velocity_function, velocity_definition.moving_time: moving_time})
dsolved_movement = dsolve(constant_accelerated_movement_definition, velocity_definition.movement_function(moving_time))
constant_accelerated_movement_function = dsolved_movement.rhs

#TODO: prove that constant_accelerated_movement_function equals to law.rhs, given C1 = initial_velocity, C2 = 0

def print():
    return pretty(law, use_unicode=False)

'''
@validate_input(charge_start_=units.charge, charge_end_=units.charge, time_=units.time)
@validate_output(units.current)
def calculate_distance(charge_start_: Quantity, charge_end_: Quantity, time_: Quantity) -> Quantity:
    charge_function_ = time * (charge_end_ - charge_start_) / time_
    applied_definition = definition.subs(charge_function(time), charge_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr, 'current')
'''