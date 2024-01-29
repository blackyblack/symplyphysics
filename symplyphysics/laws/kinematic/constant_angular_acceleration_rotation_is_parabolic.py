from sympy import Eq, solve, dsolve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    angular_acceleration_is_angular_velocity_derivative as angular_acceleration_def,
    angular_velocity_is_angle_derivative as angular_velocity_def,
)

# Description
## A body is rotating about a fixed axis with constant angular acceleration. Its angular
## displacement from initial position is a quadratic function of time and depends on
## its initial angular velocity and angular acceleration.

# Law: theta = w0 * t + alpha * t**2 / 2
## theta - angular displacement from initial position
## w0 - initial angular velocity
## alpha - constant angular acceleration
## t - time

## Conditions:
## Angular acceleration of the body is constant.

angular_displacement = Symbol("angular_displacement", angle_type)
initial_angular_velocity = Symbol("initial_angular_velocity", angle_type / units.time)
angular_acceleration = Symbol("angular_acceleration", angle_type / units.time**2)
time = Symbol("time", units.time)

law = Eq(
    angular_displacement,
    initial_angular_velocity * time + angular_acceleration * time**2 / 2,
)

def pprint(x):
    print(print_expression(x))
    assert 0

# Derive law from definitions of angular velocity and acceleration

angular_velocity_formula = dsolve(
    angular_acceleration_def.definition.subs(angular_acceleration_def.time, time),
    angular_acceleration_def.angular_velocity(time),
).rhs.subs(
    angular_acceleration_def.angular_acceleration(time), angular_acceleration,
).doit()

C1 = solve(Eq(initial_angular_velocity, angular_velocity_formula.subs(time, 0)), "C1")[0]

angular_velocity_expr = angular_velocity_formula.subs("C1", C1)

angular_displacement_formula = dsolve(
    angular_velocity_def.definition.subs(angular_velocity_def.time, time),
    angular_velocity_def.angle_function(time),
).rhs.subs(
    angular_velocity_def.angular_velocity(time), angular_velocity_expr,
).doit()

C1 = solve(Eq(0, angular_displacement_formula.subs(time, 0)), "C1")[0]

angular_displacement_expr = angular_displacement_formula.subs("C1", C1)

assert expr_equals(angular_displacement_expr, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    initial_angular_velocity_=initial_angular_velocity,
    angular_acceleration_=angular_acceleration,
    time_=time,
)
@validate_output(angular_displacement)
def calculate_angular_displacement(
    initial_angular_velocity_: Quantity,
    angular_acceleration_: Quantity,
    time_: Quantity,
) -> Quantity | float:
    result = law.rhs.subs({
        initial_angular_velocity: initial_angular_velocity_,
        angular_acceleration: angular_acceleration_,
        time: time_,
    })
    return Quantity(result)
