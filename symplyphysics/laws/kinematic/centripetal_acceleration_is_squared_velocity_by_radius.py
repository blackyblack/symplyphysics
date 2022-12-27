from symplyphysics import (
    symbols, Eq, pretty, Quantity, units, solve, simplify,
    validate_input, validate_output, expr_to_quantity, Function
)
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector
from symplyphysics.definitions import period_from_circular_frequency as circular_frequency_definition
from symplyphysics import pi, sin, cos
from sympy import diff, sqrt, atan

# Description
## When the object moves not straight (but along some curve), acceleration not only changes the quantity of velocity, but also the velocity direction.
## Velocity vector in any moment of observation is tangental to trajectory.
## In every moment of observation acceleration vector might be represented as sum of two vectors. One of them is tangental to velocity and the other one is normal and directed to center of curve.
## The tangental acceleration affects the velocity vector length and it is known as tangental acceleration, 
## and normal acceleration afects the velocity direction and also known as centripetal (that's because it's directed to center of curve).
## Law: an(t) = V(t)**2 / R
## Where:
## t is the moment of observation
## an is momental centripetal (normal) acceleration,
## V is momental linear velocity,
## R is curve radius in this point of trajectory.

centripetal_acceleration = symbols("centripetal_acceleration")
linear_velocity = symbols("linear_velocity")
curve_radius = symbols("curve_radius")

definition = Eq(centripetal_acceleration, linear_velocity**2 / curve_radius)

# Proof
## Let's assume we are having movement in 2-D space.
## Object position is described with it's radius-vector R - the vector from zero coordinates to the object and with angle alpha between X-axis and this radius-vector.
## Let's also assume we are having object is on the X-axis at the start of observation (alpha == 0, V(0) == linear_velocity), so movental velocity vector is upwards.

time = symbols("time")
alpha = symbols("alpha", cls = Function)

x_coordinate = projector.law.rhs.subs({projector.vector_length:curve_radius, projector.vector_angle: alpha(time)})
y_coordinate = projector.law.rhs.subs({projector.vector_length:curve_radius, projector.vector_angle: pi/2 - alpha(time)})

## Velocity projections are derivatives of respective coordinates.
x_velocity = diff(x_coordinate, time)
y_velocity = diff(y_coordinate, time)
velocity_vector_length = sqrt(x_velocity**2 + y_velocity**2)
## 
velocity_vector_angle = atan(y_velocity / x_velocity)

## Acceleration projections are derivatives of respective velocities.
x_acceleration = diff(x_velocity, time)
y_acceleration = diff(y_velocity, time)
full_acceleration_vector_length = sqrt(x_acceleration**2 + y_acceleration**2)
full_acceleration_vector_angle = atan(y_acceleration / x_acceleration)

angle_between_velocity_and_acceleration = full_acceleration_vector_angle - velocity_vector_angle
normal_acceleration_vector_length = simplify(full_acceleration_vector_length * sin(angle_between_velocity_and_acceleration))
print(f"{normal_acceleration_vector_length}")
## sqrt(curve_radius**2*(Derivative(alpha(time), time)**4 + Derivative(alpha(time), (time, 2))**2))*sin(atan((sin(alpha(time))*Derivative(alpha(time), time)**2 - cos(alpha(time))*Derivative(alpha(time), (time, 2)))/(sin(alpha(time))*Derivative(alpha(time), (time, 2)) + cos(alpha(time))*Derivative(alpha(time), time)**2)) + atan(1/tan(alpha(time))))
##TODO it needs little manual simplification and we'll get an = V**2 / R.

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
