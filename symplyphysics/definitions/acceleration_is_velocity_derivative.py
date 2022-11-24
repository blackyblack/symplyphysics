from symplyphysics import (
    symbols, Function, Derivative, Eq, pretty, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Acceleration definition: A = dv/dt

time = symbols('time')
acceleration, velocity_function = symbols('acceleration velocity', cls = Function)
definition = Eq(acceleration(time), Derivative(velocity_function(time), time))

definition_dimension_SI = units.meter / units.second**2

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(velocity_start_=units.velocity, velocity_end_=units.velocity, time_=units.time)
@validate_output(units.acceleration)
def calculate_linear_acceleration(velocity_start_: Quantity, velocity_end_: Quantity, time_: Quantity) -> Quantity:
    velocity_function_ = time * (velocity_end_ - velocity_start_) / time_
    applied_definition = definition.subs(velocity_function(time), velocity_function_)
    # calculate acceleration
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr, 'acceleration')
