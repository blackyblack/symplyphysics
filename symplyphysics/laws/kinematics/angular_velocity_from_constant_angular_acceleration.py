from sympy import Eq, dsolve, solve
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
    angular_acceleration_is_angular_speed_derivative as angular_acceleration_def,)

# Description
## If a body is rotating about a fixed axis with constant angular acceleration, its angular
## velocity is a linear function of time.

# Law: w = w0 + alpha * t
## w - angular velocity at time t
## w0 - initial angular velocity
## alpha - constant angular acceleration
## t - time

# Conditions
## Angular acceleration should be constant.

angular_velocity = Symbol("angular_velocity", angle_type / units.time)
initial_angular_velocity = Symbol("initial_angular_velocity", angle_type / units.time)
angular_acceleration = Symbol("angular_acceleration", angle_type / units.time**2)
time = Symbol("time", units.time)

law = Eq(angular_velocity, initial_angular_velocity + angular_acceleration * time)

# Derive this law from definition of angular acceleration

angular_velocity_formula = dsolve(
    angular_acceleration_def.definition.subs(angular_acceleration_def.time, time),
    angular_acceleration_def.angular_speed(time),
).rhs.subs(
    angular_acceleration_def.angular_acceleration(time),
    angular_acceleration,
).doit()

angular_velocity_derived = solve([
    Eq(initial_angular_velocity, angular_velocity_formula.subs(time, 0)),
    Eq(angular_velocity, angular_velocity_formula)
], ("C1", angular_velocity),
    dict=True)[0][angular_velocity]

assert expr_equals(angular_velocity_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    initial_angular_velocity_=initial_angular_velocity,
    angular_acceleration_=angular_acceleration,
    time_=time,
)
@validate_output(angular_velocity)
def calculate_angular_velocity(
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
