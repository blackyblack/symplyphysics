from symplyphysics import (
    symbols, Eq, pretty, Quantity, units, solve, simplify,
    validate_input, validate_output, expr_to_quantity, Function, Derivative,
    pi, sin, cos, diff, sqrt, atan
)
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector

# Description
## When the object moves not straight but along some curve, acceleration not only changes the magnitude of velocity, but also the velocity direction.
## Velocity vector in any moment of observation is tangential to trajectory.
## In every moment of observation acceleration vector might be represented as sum of two vectors. One of them is tangential to velocity and the other one is normal and directed to the center of curve.
## The tangential acceleration affects the velocity vector length and it is known as tangential acceleration. 
## Normal acceleration afects the velocity direction and also known as centripetal (that's because it's directed to the center of curve).

# Law: an(t) = V(t)**2 / R
## Where:
## t is the moment of observation
## an is momental centripetal (normal) acceleration
## V is momental linear velocity
## R is curve radius in this point of trajectory.

centripetal_acceleration = symbols("centripetal_acceleration")
linear_velocity = symbols("linear_velocity")
curve_radius = symbols("curve_radius")

law = Eq(centripetal_acceleration, linear_velocity**2 / curve_radius)

# Derive the same law from acceleration and velocity definitions

## Let's assume we are having movement in 2-D space.
## Object position is described with it's radius-vector R - the vector from zero coordinates to the object and with angle 'alpha' between X-axis and this radius-vector.

time = symbols("time")
alpha = symbols("alpha", cls = Function, positive=True)

curve_radius_horisontal = projector.law.rhs.subs({projector.vector_length: curve_radius, projector.vector_angle: alpha(time)})
curve_radius_vertical = projector.law.rhs.subs({projector.vector_length: curve_radius, projector.vector_angle: pi/2 - alpha(time)})

## Velocity projections are derivatives of respective coordinates.

#TODO: use velocity definition as soon as it is added to symplyphysics
velocity_horisontal = diff(curve_radius_horisontal, time)
velocity_vertical = diff(curve_radius_vertical, time)
velocity_vector_length = sqrt(velocity_horisontal**2 + velocity_vertical**2)
#TODO: explain why
velocity_vector_angle = atan(velocity_vertical / velocity_horisontal)

## velocity_vector_angle = -atan(cos(alpha(time))/sin(alpha(time)))
## velocity_vector_angle = -atan(sin(pi/2 - alpha(time)) / cos(pi/2 - alpha(time))) = -atan(tan(pi/2 - alpha(time)))
## velocity_vector_angle = alpha(time) - pi/2, for alpha(time) in range (0, pi)
velocity_vector_angle = alpha(time) - pi/2

## Acceleration projections are derivatives of respective velocities.

velocity_vector = [velocity_horisontal, velocity_vertical]
#TODO: derive from projector
radial_unit_vector = [cos(alpha(time)), sin(alpha(time))]
tangential_unit_vector = [-sin(alpha(time)), cos(alpha(time))]

## Tangential unit vector is radial unit vector derivative
#TODO: use Derivative operator and Vector class
assert diff(radial_unit_vector[0], alpha(time)) == tangential_unit_vector[0]
assert diff(radial_unit_vector[1], alpha(time)) == tangential_unit_vector[1]

## This is Dot product of radial vector and velocity vector. Radial vector is orthogonal to velocity hence vector
## multiplication result should be zero.
#TODO: use Dot operator and Vector class
assert velocity_vector[0] * radial_unit_vector[0] + velocity_vector[1] * radial_unit_vector[1] == 0

#TODO: use acceleration definition
acceleration_horisontal = diff(velocity_horisontal, time)
acceleration_vertical = diff(velocity_vertical, time)
acceleration_vector = [acceleration_horisontal, acceleration_vertical]

tangential_acceleration_magnitude = curve_radius * Derivative(alpha(time), (time, 2))
radial_acceleration_magnitude = -curve_radius * Derivative(alpha(time), time)**2
#TODO: use scalar multiplication and Vector class
tangential_acceleration = [tangential_acceleration_magnitude * tangential_unit_vector[0], tangential_acceleration_magnitude * tangential_unit_vector[1]]
radial_acceleration = [radial_acceleration_magnitude * radial_unit_vector[0], radial_acceleration_magnitude * radial_unit_vector[1]]

#TODO: use vector addition and Vector class
assert simplify(acceleration_vector[0] - (tangential_acceleration[0] + radial_acceleration[0])) == 0
assert simplify(acceleration_vector[1] - (tangential_acceleration[1] + radial_acceleration[1])) == 0

## Here we've proven that tangential_acceleration + radial_acceleration equals to acceleration_vector. It means, we've
## changed basis of acceleration_vector to tangential and radial vectors instead of cartesian coordinates.
## Same result could be achieved by rotating coordinate system to 'velocity_vector_angle'

## We are not interested in tangential_acceleration as we are looking for centripetal acceleration which is 'radial_acceleration'
## in our proof.

#TODO: add law with linear_velocity == Derivative(alpha(time), time) * curve_radius
law_acceleration = law.rhs.subs(linear_velocity, Derivative(alpha(time), time) * curve_radius)

## radial_acceleration_magnitude has minus sign. It means it is directed towards the center of the curve. The centripetal
## acceleration law is not defined in vector terms so we should only compare acceleration magnitudes (absolute values).
## And for some reason, sympy does not allow to compare Abs with non-Abs values so we apply abs() to both sides.
assert simplify(abs(radial_acceleration_magnitude) - abs(law_acceleration)) == 0


def print():
    return pretty(law, use_unicode=False)

@validate_input(linear_velocity_=units.velocity, curve_radius_=units.length)
@validate_output(units.acceleration)
def calculate_acceleration(linear_velocity_: Quantity, curve_radius_: Quantity) -> Quantity:        
    solved = solve(law, centripetal_acceleration, dict=True)[0][centripetal_acceleration]
    result_expr = solved.subs({
        linear_velocity:linear_velocity_,
        curve_radius:curve_radius_
    })
    return expr_to_quantity(result_expr, "centripetal_acceleration")
