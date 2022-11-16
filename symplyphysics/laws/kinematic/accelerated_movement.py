from symplyphysics import (
    symbols, Function, Derivative, Eq, pretty, Quantity, units,
    validate_input, validate_output, expr_to_quantity
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
constant_acceleration = symbols('constant_acceleration')
initial_velocity = symbols('initial_velocity')
theoretical_distance_function, distance_function, velocity_function, acceleration_function = symbols('theoretical_distance_function distance_function velocity_function acceleration_function', cls = Function)

law = Eq(theoretical_distance_function(moving_time), initial_velocity * moving_time + constant_acceleration * moving_time * moving_time / 2)

constant_acceleration_definition = acceleration_definition.definition.subs({acceleration_definition.acceleration: constant_acceleration, acceleration_definition.time: moving_time})
dsolved_velocity = constant_acceleration_definition.doit()
constant_accelerated_velocity_function = dsolved_velocity.rhs

constant_accelerated_movement_definition = velocity_definition.definition.subs(velocity_function(moving_time), constant_accelerated_velocity_function)
dsolved_movement = constant_accelerated_movement_definition.doit()
constant_accelerated_movement_function = dsolved_movement.rhs

definition_dimension_SI = units.meter

def print1():
    return pretty(constant_accelerated_velocity_function, use_unicode=False)

def print():
    return pretty(constant_accelerated_movement_function, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)



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