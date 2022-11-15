from symplyphysics import (
    symbols, Function, Derivative, Eq, pretty, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Accelerated movement is the kind of movement when object has constant acceleration (e.g with the constant force applied to object).
## For example, free fall is accelerated movement. The acceleration in each moment of time is g.
## According to velocity and acceleration definitions, velocity is movement derivative, acceleration is velocity derivative.

moving_time = symbols('moving_time')
constant_acceleration = symbols('constant_acceleration')
distance_function, velocity_function, acceleration_function = symbols('distance_function velocity_function acceleration_function', cls = Function)

velocity_definition = Eq(velocity_function(moving_time), Derivative(distance_function(moving_time), moving_time))
acceleration_definition = Eq(acceleration_function(moving_time), Derivative(velocity_function(moving_time), moving_time))

constant_acceleration_definition = acceleration_definition.subs(acceleration_function(moving_time), constant_acceleration)
dsolved_velocity = constant_acceleration_definition.doit()
constant_accelerated_velocity_function = dsolved_velocity.rhs

constant_accelerated_movement_definition = velocity_definition.subs(velocity_function(moving_time), constant_accelerated_velocity_function)
dsolved_movement = constant_accelerated_movement_definition.doit()
constant_accelerated_movement_function = dsolved_movement.rhs

definition_dimension_SI = units.meter

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