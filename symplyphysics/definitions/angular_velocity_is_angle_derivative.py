from symplyphysics import (
    symbols, Function, Derivative, Eq, pretty, Quantity, units, pi, SI,
    validate_input, validate_output, expr_to_quantity
)

from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

# Description
## The movement along circle might be easily represented in polar coordinates with the pole in the center of the circle.
## Object position is described by radius-vector from pole to object. Length of this vector is radius of circle and
## angle of this vector and it is also known as movement phase.

# Definition: ω(t) = φ(t)/dt, where
## ω(t) is angular velocity function of time
## φ(t) is angle function of time

time = symbols('time')
angular_velocity_function, angle_function = symbols('angular_velocity_function angle_function', cls = Function)

definition = Eq(angular_velocity_function(time), Derivative(angle_function(time), time))

definition_dimension_SI = units.radian / units.second

#def print():
#    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(phase_start_= angle_type, phase_end_= angle_type, moving_time_=units.time)
@validate_output(angle_type / units.time)
def calculate_velocity(phase_start_: Quantity, phase_end_: Quantity, moving_time_: Quantity) -> Quantity:
    phase_function_ = time * (phase_end_ - phase_start_) / moving_time_
    print(f"{phase_function_}")
    applied_definition = definition.subs(angle_function(time), phase_function_)
    print(f"{applied_definition}")
    dsolved = applied_definition.doit()
    print(f"{dsolved}")
    result_expr = dsolved.rhs
    print(f"{result_expr}")
    result_qty = expr_to_quantity(result_expr, 'angular_velocity')
    print(f"{result_qty}")
    return result_qty


a0 = units.Quantity('a0')
SI.set_quantity_dimension(a0, angle_type)
SI.set_quantity_scale_factor(a0, 0 * units.radian)
a1 = units.Quantity('a1')
SI.set_quantity_dimension(a1, angle_type)
SI.set_quantity_scale_factor(a1, pi * units.radian)    
t = units.Quantity('t')
SI.set_quantity_dimension(t, units.time)
SI.set_quantity_scale_factor(t, 5 * units.second)

calculate_velocity(a0, a1, t)
