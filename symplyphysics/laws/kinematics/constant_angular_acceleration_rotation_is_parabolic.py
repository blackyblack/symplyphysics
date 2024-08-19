from sympy import Eq, solve, dsolve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    angular_acceleration_is_angular_speed_derivative as angular_acceleration_def,
    angular_speed_is_angular_distance_derivative as angular_velocity_def,
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
## - Angular acceleration of the body is constant.
## - Initial displacement is zero.

angular_displacement = Symbol("angular_displacement", angle_type)
initial_angular_velocity = Symbol("initial_angular_velocity", angle_type / units.time)
angular_acceleration = Symbol("angular_acceleration", angle_type / units.time**2)
time = Symbol("time", units.time)

law = Eq(
    angular_displacement,
    initial_angular_velocity * time + angular_acceleration * time**2 / 2,
)

# Derive law from definitions of angular velocity and acceleration

_angular_velocity_formula = dsolve(
    angular_acceleration_def.definition.subs(angular_acceleration_def.time, time),
    angular_acceleration_def.angular_speed(time),
).rhs.subs(
    angular_acceleration_def.angular_acceleration(time),
    angular_acceleration,
).doit()

_angular_velocity = Symbol("_angular_velocity", angle_type / units.time)
_angular_velocity_derived = solve(
    [
    Eq(initial_angular_velocity, _angular_velocity_formula.subs(time, 0)),
    Eq(_angular_velocity, _angular_velocity_formula)
    ],
    ("C1", _angular_velocity),
    dict=True,
)[0][_angular_velocity]

_angular_displacement_formula = dsolve(
    angular_velocity_def.definition.subs(angular_velocity_def.time, time),
    angular_velocity_def.angular_distance(time),
).rhs.subs(
    angular_velocity_def.angular_speed(time),
    _angular_velocity_derived,
).doit()

_angular_displacement_derived = solve(
    [
    # initial angular displacement is 0 by condition
    Eq(0, _angular_displacement_formula.subs(time, 0)),
    Eq(angular_displacement, _angular_displacement_formula)
    ],
    ("C1", angular_displacement),
    dict=True,
)[0][angular_displacement]

assert expr_equals(_angular_displacement_derived, law.rhs)


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
) -> Quantity:
    result = law.rhs.subs({
        initial_angular_velocity: initial_angular_velocity_,
        angular_acceleration: angular_acceleration_,
        time: time_,
    })
    return Quantity(result)
