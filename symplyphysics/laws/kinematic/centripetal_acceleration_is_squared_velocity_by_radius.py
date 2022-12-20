from symplyphysics import (
    symbols, Eq, pretty, Quantity, units, solve, simplify,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector
from symplyphysics.definitions import period_from_circular_frequency as circular_frequency_definition
from symplyphysics import pi, sin, cos
from sympy import diff, sqrt, atan

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
## Object position is described with it's radius-vector R - the vector from zero coordinates to the object and with angle between X-axis and this radius-vector.
## Tgis angle is called phase. 
## The angle between Y-axis and radius-vector is (pi/2 - phase).
## Let center of circle have zero coordinates. 
## Let the alpha be 0 at the start of observation.

phase = symbols("phase")
observation_time = symbols("_observation_time")
period = symbols("period")

## Phase is the function of time. Circular velocity is the derivative of phase by time. If the spinning is uniform (harmonic) with period of T, circular velocity is constant.
circular_frequency = solve(circular_frequency_definition.definition, circular_frequency_definition.circular_frequency, dict = True)[0][circular_frequency_definition.circular_frequency]
phase = circular_frequency * observation_time

x_coordinate = projector.law.rhs.subs({projector.vector_length:curve_radius, projector.vector_angle: phase})
y_coordinate = projector.law.rhs.subs({projector.vector_length:curve_radius, projector.vector_angle: pi/2 - phase})
## BTW x and y at any moment of time satisfy equation x**2 + y**2 == R**2, but i don't know if it is gonna be useful here or not.

## Velocity projections are derivatives of respective coordinates.
x_velocity = diff(x_coordinate, observation_time)
y_velocity = diff(y_coordinate, observation_time)
velocity_vector_length = simplify(sqrt(x_velocity**2 + y_velocity**2))
velocity_vector_angle = atan(y_velocity / x_velocity)
## Velocity vector length is independent from alpha:
assert(velocity_vector_length == 2*pi*sqrt(curve_radius**2/period**2))

## Linear velocity vector angle is -(pi/2 - phase)
assert(velocity_vector_angle == -atan(cos(2*pi*observation_time/period)/sin(2*pi*observation_time/period)))

##TODO with properly defined circular velocity it is all too difficult for understanding. let's define it first.


## Acceleration projections are derivatives of respective velocities.
x_acceleration = diff(x_velocity, observation_time)
y_acceleration = diff(y_velocity, observation_time)

acceleration_vector_length = simplify(sqrt(x_acceleration**2 + y_acceleration**2))
print(f"{acceleration_vector_length}")
acceleration_angle = atan(y_acceleration / x_acceleration)

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
