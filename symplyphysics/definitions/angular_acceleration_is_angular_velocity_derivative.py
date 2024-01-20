from sympy import (Eq, Derivative)
from symplyphysics import (angle_type, units, Quantity, Function, Symbol, print_expression,
    validate_input, validate_output)

# Description
## Angular acceleration is a kinematic quantity that characterizes the change in angular velocity over time.
## The word "kinematic" means that motion is considered without taking
## into account the action of forces on the body, regardless of them.

# Definition: epsilon(t) = ω(t)/dt
# Where:
## ω(t) is angular velocity function of time
## epsilon(t) is anglular acceleration of time

time = Symbol("time", units.time)
angular_velocity = Function("angular_velocity", angle_type / units.time)
angular_acceleration = Function("angular_acceleration", angle_type / (units.time ** 2))

definition = Eq(angular_acceleration(time), Derivative(angular_velocity(time), time))

definition_units_SI = units.radian / (units.second ** 2)


def print_law() -> str:
    return print_expression(definition)


@validate_input(angular_velocity_start_=angular_velocity, angular_velocity_end_=angular_velocity, moving_time_=time)
@validate_output(angular_acceleration)
def calculate_angular_acceleration(angular_velocity_start_: Quantity,
    angular_velocity_end_: Quantity, moving_time_: Quantity) -> Quantity:
    #HACK: SymPy angles are always in radians
    angular_velocity_function = time * (angular_velocity_end_ - angular_velocity_start_) / moving_time_
    applied_definition = definition.subs(angular_velocity(time), angular_velocity_function)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
