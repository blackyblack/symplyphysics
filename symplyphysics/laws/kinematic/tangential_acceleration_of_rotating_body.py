from sympy import Eq, solve, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    angle_type
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic import linear_velocity_from_angular_velocity_and_radius as linear_velocity_law

# Description
## The tangential acceleration of a rotating body represents the change in magnitude
## of the velocity vector, and is tangent to the path of the body. It is proportional
## to the angular acceleration of the body and its rotation radius.

# Law: a_t = alpha * r
## a_t - tangential acceleration
## alpha - angular acceleration
## r - rotation radius

tangential_acceleration = Symbol("tangential_acceleration", units.acceleration)
angular_acceleration = Symbol("angular_acceleration", angle_type / units.time**2)
rotation_radius = Symbol("rotation_radius", units.length)

law = Eq(tangential_acceleration, angular_acceleration * rotation_radius)


# Derive the law from the law of linear velocity in case of a rotating body
# Condition: radius of rotation remains constant.

linear_velocity = Function("linear_velocity", units.velocity)
angular_velocity = Function("angular_velocity", angle_type / units.time)
time = Symbol("time", units.time)

linear_velocity_law_sub = linear_velocity_law.law.subs({
    linear_velocity_law.linear_velocity: linear_velocity(time),
    linear_velocity_law.angular_velocity: angular_velocity(time),
    linear_velocity_law.curve_radius: rotation_radius,
})

# Differentiate both sides of the equation w.r.t. time
diff_linear_velocity_law = Eq(
    linear_velocity_law_sub.lhs.diff(time),
    linear_velocity_law_sub.rhs.diff(time),
)

# alpha = d(omega)/dt
diff_linear_velocity_law_sub = diff_linear_velocity_law.subs(
    Derivative(angular_velocity(time), time), 
    angular_acceleration
)

# a_t = dv/dt
tangential_acceleration_derived = solve(
    diff_linear_velocity_law_sub,
    Derivative(linear_velocity(time), time)
)[0]

tangential_acceleration_from_law = solve(law, tangential_acceleration)[0]

assert expr_equals(
    tangential_acceleration_derived,
    tangential_acceleration_from_law,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    angular_acceleration_=angular_acceleration,
    rotation_radius_=rotation_radius,
)
@validate_output(tangential_acceleration)
def calculate_tangential_acceleration(
    angular_acceleration_: Quantity, 
    rotation_radius_: Quantity
) -> Quantity:
    result_expr = solve(law, tangential_acceleration)[0]
    result = result_expr.subs({
        angular_acceleration: angular_acceleration_,
        rotation_radius: rotation_radius_,
    })
    return Quantity(result)
