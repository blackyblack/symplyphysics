from symplyphysics import (
    symbols, Eq, pretty, Quantity, units, solve, SI,
    validate_input, validate_output, expr_to_quantity
)
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

# Description
## Angular velocity is the rate of change of the angular position of a rotating body. We can define the angular velocity of a particle as the rate
## at which the particle rotates around a centre point i.e., the time rate of change of its angular displacement relative to the origin. Linear
## velocity is the measure of the rate of change of displacement with respect to time when the object moves along a straight path.

# Law: V = ω * R
## Where:
## V is momental linear velocity
## ω is angular velocity
## R is curve radius in this point of trajectory.

linear_velocity, angular_velocity = symbols("linear_velocity angular_velocity")
curve_radius = symbols("curve_radius")

law = Eq(linear_velocity, angular_velocity * curve_radius)


def print():
    return pretty(law, use_unicode=False)

@validate_input(angular_velocity_=(angle_type / units.time), curve_radius_=units.length)
@validate_output(units.velocity)
def calculate_linear_velocity(angular_velocity_: Quantity, curve_radius_: Quantity) -> Quantity:        
    solved = solve(law, linear_velocity, dict=True)[0][linear_velocity]
    #HACK: sympy angles are always in radians and angular velocity cannot be properly converted to velocity
    SI.set_quantity_dimension(angular_velocity_, 1 / units.time)
    result_expr = solved.subs({angular_velocity: angular_velocity_, curve_radius: curve_radius_})
    return expr_to_quantity(result_expr, "linear_velocity_out")
