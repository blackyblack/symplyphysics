from sympy import (Eq, solve, sin, cos, Derivative, pi)
from symplyphysics import (units, Quantity, Symbol, Function, print_expression, angle_type,
    CoordinateSystem, Vector, validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals, expr_equals_abs
from symplyphysics.core.vectors.arithmetics import dot_vectors
from symplyphysics.definitions import velocity_is_movement_derivative as velocity_def
from symplyphysics.definitions import angular_velocity_is_angle_derivative as angular_velocity_def
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration_def
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector
from symplyphysics.laws.kinematic import linear_velocity_from_angular_velocity_and_radius as linear_velocity_law

# Description
## When the object moves not straight but along some curve, acceleration not only changes the magnitude of velocity, but also the velocity direction.
## Velocity vector in any moment of observation is tangential to trajectory.
## In every moment of observation acceleration vector might be represented as sum of two vectors. One of them is tangential to velocity and the other one is normal and directed to the center of curve.
## The tangential acceleration affects the velocity vector length and it is known as tangential acceleration.
## Normal (radial) acceleration afects the velocity direction and also known as centripetal (that's because it's directed to the center of curve).
## 3-D and higher dimensional motion also have centripetal acceleration but a generalized form of the law should be applied.

# Law: an(t) = V(t)**2 / R
## Where:
## t is the moment of observation
## an is momental centripetal (normal) acceleration
## V is momental linear velocity
## R is curve radius in this point of trajectory.

# Conditions
## - Curve is smooth and continuous. Does not have 0 radius.
## - Motion is in 2-D space.

centripetal_acceleration = Symbol("centripetal_acceleration", units.acceleration)
linear_velocity = Symbol("linear_velocity", units.velocity)
curve_radius = Symbol("curve_radius", units.length)

law = Eq(centripetal_acceleration, linear_velocity**2 / curve_radius)

# Derive the same law from acceleration and velocity definitions

## Let's assume we are having movement in 2-D space.
## Object position is described with it's radius-vector R - the vector from zero coordinates to the object and with angle 'alpha' between X-axis and this radius-vector.

time = Symbol("time", units.time)
alpha = Function("alpha", angle_type, positive=True)
cartesian_coordinates = CoordinateSystem()

curve_radius_horisontal = projector.law.rhs.subs({
    projector.vector_length: curve_radius,
    projector.vector_angle: alpha(time)
})
curve_radius_vertical = projector.law.rhs.subs({
    projector.vector_length: curve_radius,
    projector.vector_angle: pi / 2 - alpha(time)
})

## Velocity projections are derivatives of respective coordinates.

#NOTE: replace 'moving_time' first as Derivative can have difficulties when processing both substitutions at once
velocity_horisontal = velocity_def.definition.rhs.subs(velocity_def.moving_time,
    time).subs(velocity_def.movement(time), curve_radius_horisontal).doit()
velocity_vertical = velocity_def.definition.rhs.subs(velocity_def.moving_time,
    time).subs(velocity_def.movement(time), curve_radius_vertical).doit()
velocity_vector = Vector([velocity_horisontal, velocity_vertical], cartesian_coordinates)

## These unit vectors should not necessary be derived. We can choose them at will and prove that
## they are orthogonal to each other and radial_unit_vector is orthogonal to 'velocity_vector'.
## One can also show that 'tangential_unit_vector' is 'radial_unit_vector' derivative.
radial_unit_vector = Vector([cos(alpha(time)), sin(alpha(time))], cartesian_coordinates)
tangential_unit_vector = Vector([-sin(alpha(time)), cos(alpha(time))], cartesian_coordinates)

## This is Dot product of radial vector and velocity vector. Radial vector is orthogonal to velocity hence vector
## multiplication result should be zero.
assert expr_equals(dot_vectors(radial_unit_vector, velocity_vector), 0)
## Radial vector is orthogonal to tangential vector hence tangential vector should be parallel to velocity vector.
assert expr_equals(dot_vectors(tangential_unit_vector, radial_unit_vector), 0)

## Use acceleration definition to calculate 'acceleration_vector'
acceleration_horisontal = acceleration_def.definition.rhs.subs(acceleration_def.time, time)
acceleration_horisontal = acceleration_horisontal.subs(acceleration_def.velocity(time),
    velocity_horisontal).doit()
acceleration_vertical = acceleration_def.definition.rhs.subs(acceleration_def.time, time)
acceleration_vertical = acceleration_vertical.subs(acceleration_def.velocity(time),
    velocity_vertical).doit()
acceleration_vector = Vector([acceleration_horisontal, acceleration_vertical],
    cartesian_coordinates)

## Prove that 'acceleration_vector' has tangential and radial parts.

tangential_acceleration_magnitude = curve_radius * Derivative(alpha(time), (time, 2))
radial_acceleration_magnitude = -curve_radius * Derivative(alpha(time), time)**2

## Use Dot product to find tangential and radial components of acceleration. Confirm they are
## equal to expected value: tangential_acceleration_magnitude, radial_acceleration_magnitude
tangential_acceleration_component = dot_vectors(acceleration_vector, tangential_unit_vector)
radial_acceleration_component = dot_vectors(acceleration_vector, radial_unit_vector)
assert expr_equals(tangential_acceleration_component, tangential_acceleration_magnitude)
assert expr_equals(radial_acceleration_component, radial_acceleration_magnitude)

## Here we've proven that tangential_acceleration + radial_acceleration equals to acceleration_vector. It means, we've
## changed basis of acceleration_vector to tangential and radial vectors instead of cartesian coordinates.
## Same result could be achieved by rotating coordinate system by velocity vector angle.

## We are not interested in tangential_acceleration as we are looking for centripetal acceleration which is 'radial_acceleration'
## in our proof.

angular_velocity_applied = angular_velocity_def.definition.rhs.subs(angular_velocity_def.time, time)
angular_velocity_applied = angular_velocity_applied.subs(angular_velocity_def.angle_function(time),
    alpha(time))
linear_velocity_applied = linear_velocity_law.law.rhs.subs({
    linear_velocity_law.angular_velocity: angular_velocity_applied,
    linear_velocity_law.curve_radius: curve_radius
})
law_acceleration = law.rhs.subs(linear_velocity, linear_velocity_applied)

## radial_acceleration_magnitude has minus sign. It means it is directed towards the center of the curve. The centripetal
## acceleration law is not defined in vector terms so we should only compare acceleration magnitudes (absolute values).
assert expr_equals_abs(radial_acceleration_magnitude, law_acceleration)


def print_law() -> str:
    return print_expression(law)


@validate_input(linear_velocity_=linear_velocity, curve_radius_=curve_radius)
@validate_output(centripetal_acceleration)
def calculate_acceleration(linear_velocity_: Quantity, curve_radius_: Quantity) -> Quantity:
    solved = solve(law, centripetal_acceleration, dict=True)[0][centripetal_acceleration]
    result_expr = solved.subs({linear_velocity: linear_velocity_, curve_radius: curve_radius_})
    return Quantity(result_expr)
