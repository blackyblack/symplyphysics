from symplyphysics import (
    symbols, Eq, pretty, Quantity, units, solve, diff, sqrt, simplify, arctan,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector
from symplyphysics import pi


# Description
## Centripetal acceleration is type of acceleration which is perpendicular to velocity vector and directed to the curve center.
## Centripetal acceleration only changes the velocity vector direction and not the velocity vector length.

## Law: an = V**2 / R
## Where:
## an is centripetal acceleration (aka normal acceleration),
## V is linear velocity,
## R is curve radius.

centripetal_acceleration = symbols("centripetal_acceleration")
linear_velocity = symbols("linear_velocity")
curve_radius = symbols("curve_radius")

definition = Eq(centripetal_acceleration, linear_velocity**2 / curve_radius)

## Circular movement might be described with 2-dimensional coordinate system as any circle is 2-D.
## Object position is described with it's radius-vector R - the vector from zero coordinates to the object and with angle alpha between X-axis and this radius-vector. 
## The angle between Y-axis and radius-vector is (pi/2 - alpha).
## Let center of circle have zero coordinates. 
## Let the alpha be 0 at the start of observation.

alpha = symbols("alpha")

x_coordinate = projector.law.rhs.subs({projector.vector_length:curve_radius, projector.vector_angle: alpha})
y_coordinate = projector.law.rhs.subs({projector.vector_length:curve_radius, projector.vector_angle: pi/2 - alpha})
## BTW x and y at any moment of time satisfy equation x**2 + y**2 == R**2, but i don't know if it is useful here or not.

## Velocity projections are derivatives of respective coordinates.
x_velocity = diff(x_coordinate, alpha)
y_velocity = diff(y_coordinate, alpha)
velocity_vector_length = simplify(sqrt(x_velocity**2 + y_velocity**2))
velocity_vector_angle = arctan(y_velocity / x_velocity)
## Velocity vector length is independent from alpha

## Acceleration projections are derivatives of respective velocities.
x_acceleration = diff(x_velocity, alpha)
y_acceleration = diff(y_velocity, alpha)

total_acceleration = sqrt(x_acceleration**2 + y_acceleration**2)
acceleration_angle = arctan(y_acceleration / x_acceleration)

## We can prove angle between velocity and acceleration is pi/2
## also we can prove velocity is tangent.

def print():
    return pretty(definition, use_unicode=False)

@validate_input(linear_velocity_=units.velocity, curve_radius_=units.length)
@validate_output(units.acceleration)
def calculate_acceleration(linear_velocity_: Quantity, curve_radius_: Quantity) -> Quantity:        
    solved = solve(definition, centripetal_acceleration, dict=True)[0][centripetal_acceleration]
    result_expr = solved.subs({
        linear_velocity:linear_velocity_,
        curve_radius:curve_radius_
    })
    return expr_to_quantity(result_expr, "acceleration")
