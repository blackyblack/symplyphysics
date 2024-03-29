from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, angle_type, validate_input,
    validate_output)

# Description
## Angular velocity is the rate of change of the angular position of a rotating body. We can define the angular velocity of a particle as the rate
## at which the particle rotates around a centre point i.e., the time rate of change of its angular displacement relative to the origin. Linear
## velocity is the measure of the rate of change of displacement with respect to time when the object moves along a straight path.

# Law: V = ω * R
## Where:
## V is momental linear velocity
## ω is angular velocity
## R is curve radius in this point of trajectory.

# Note
## - This formula also applies to the speed of the center of mass of a wheel that is rolling smoothly.

linear_velocity = Symbol("linear_velocity", units.velocity)
angular_velocity = Symbol("angular_velocity", angle_type / units.time)
curve_radius = Symbol("curve_radius", units.length)

law = Eq(linear_velocity, angular_velocity * curve_radius)


def print_law() -> str:
    return print_expression(law)


@validate_input(angular_velocity_=angular_velocity, curve_radius_=curve_radius)
@validate_output(linear_velocity)
def calculate_linear_velocity(angular_velocity_: Quantity, curve_radius_: Quantity) -> Quantity:
    solved = solve(law, linear_velocity, dict=True)[0][linear_velocity]
    result_expr = solved.subs({angular_velocity: angular_velocity_, curve_radius: curve_radius_})
    return Quantity(result_expr)
