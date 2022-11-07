from symplyphysics import (
    symbols, Function, Derivative, Eq, pretty, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Velocity V(t) = dS(t)/dt
## Where V(t) is velocity function of time, S(t) is movement

time = symbols('time')
velocity_function, movement_function = symbols('velocity movement', cls = Function)
definition = Eq(velocity_function(time), Derivative(movement_function(time), time))

definition_dimension_SI = units.meter / units.second

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(position_start_=units.length, position_end_=units.length, time_=units.time)
@validate_output(units.meter / units.second)
def calculate_velocity(position_start_: Quantity, position_end_: Quantity, time_: Quantity) -> Quantity:
    movement_function_ = time * (position_end_ - position_start_) / time_
    applied_definition = definition.subs(movement_function(time), movement_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr, 'velocity')
