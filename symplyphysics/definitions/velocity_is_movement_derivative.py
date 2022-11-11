from symplyphysics import (
    symbols, Function, Derivative, Eq, pretty, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Kinematical velocity definition: V(t) = dS(t)/dt, where
## V(t) is velocity function from time
## S(t) is movement function from time

moving_time = symbols('moving_time')
velocity_function, movement_function = symbols('current charge', cls = Function)
definition = Eq(velocity_function(moving_time), Derivative(movement_function(moving_time), moving_time))

definition_dimension_SI = units.meter / units.second

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(movement_start_=units.length, movement_end_=units.length, moving_time_=units.time)
@validate_output(units.velocity)
def calculate_velocity(movement_start_: Quantity, movement_end_: Quantity, moving_time_: Quantity) -> Quantity:
    movement_function_ = moving_time * (movement_end_ - movement_start_) / moving_time_
    applied_definition = definition.subs(movement_function(moving_time), movement_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr, 'velocity')
