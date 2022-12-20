from symplyphysics import (
    symbols, Function, Derivative, Eq, pretty, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

# Description
## The movement along circle might be easily represented in polar coordinates with the pole in the center of the circle.
## Object position is described by radius-vector from pole to object. Length of this vector is radius of circle and
## angle of this vector is movement phase.

# Definition: ω(t) = φ(t)/dt, where
## ω(t) is circular velocity function of time
## φ(t) is phase function of time

observation_time = symbols('observation_time')
circular_velocity_function, phase_function = symbols('circular_velocity_function phase_function', cls = Function)
definition = Eq(circular_velocity_function(observation_time), Derivative(phase_function(observation_time), observation_time))

definition_dimension_SI = units.radian / units.second

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(phase_start_=units.angle_type, position_end_=units.angle_type, moving_time_=units.time)
@validate_output(units.angle_type / units.time)
def calculate_velocity(position_start_: Quantity, position_end_: Quantity, moving_time_: Quantity) -> Quantity:
    phase_function_ = observation_time * (position_end_ - position_start_) / moving_time_
    applied_definition = definition.subs(phase_function(observation_time), phase_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr, 'circular_velocity')
