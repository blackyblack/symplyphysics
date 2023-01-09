from symplyphysics import (
    symbols, Function, Derivative, Eq, pretty, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

# Description
## The movement along circle might be easily represented in polar coordinates with the pole in the center of the circle.
## Object position is described by radius-vector from pole to object. 
## Magnitude of this vector is radius of circle. Angle between this vector and X axis it is also known as movement phase.

# Definition: ω(t) = φ(t)/dt, where
## ω(t) is angular velocity function of time
## φ(t) is angle function of time

time = symbols('time')
angular_velocity_function, angle_function = symbols('angular_velocity_function angle_function', cls = Function)

definition = Eq(angular_velocity_function(time), Derivative(angle_function(time), time))

definition_dimension_SI = units.radian / units.second

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(angle_start_= angle_type, angle_end_= angle_type, moving_time_=units.time)
@validate_output(1 / units.time)
def calculate_angular_velocity(angle_start_: Quantity, angle_end_: Quantity, moving_time_: Quantity) -> Quantity:
    #HACK: sympy angles are always in radians
    angle_start_radians = angle_start_.scale_factor
    angle_end_radians = angle_end_.scale_factor
    angle_function_ = time * (angle_end_radians - angle_start_radians) / moving_time_
    applied_definition = definition.subs(angle_function(time), angle_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr, "angular_velocity")
