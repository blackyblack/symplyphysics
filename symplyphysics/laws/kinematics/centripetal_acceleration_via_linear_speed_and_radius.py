"""
Centripetal acceleration via linear speed and radius
====================================================

*Centripetal acceleration* is defined as the change in velocity tangential to the velocity vector.

**Links:**

#. `Physics LibreTexts, first part of equation 6.2.5 <https://phys.libretexts.org/Bookshelves/College_Physics/College_Physics_1e_(OpenStax)/06%3A_Uniform_Circular_Motion_and_Gravitation/6.02%3A_Centripetal_Acceleration>`__.
"""

from sympy import Eq, solve, sin, cos, Derivative, pi
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    CoordinateSystem,
    Vector,
    validate_input,
    validate_output,
    clone_as_function,
)
from symplyphysics.core.expr_comparisons import expr_equals, expr_equals_abs
from symplyphysics.core.vectors.arithmetics import dot_vectors
from symplyphysics.definitions import speed_is_distance_derivative as velocity_def
from symplyphysics.definitions import angular_speed_is_angular_distance_derivative as angular_velocity_def
from symplyphysics.definitions import acceleration_is_speed_derivative as acceleration_def
from symplyphysics.laws.geometry import planar_projection_is_cosine as projector
from symplyphysics.laws.kinematics import speed_via_angular_speed_and_radius as linear_velocity_law

centripetal_acceleration = clone_as_symbol(symbols.acceleration, subscript="n")
"""
Centripetal, or normal, :symbols:`acceleration`.
"""

speed = symbols.speed
"""
Linear :symbols:`speed`.
"""

radius_of_curvature = symbols.radius_of_curvature
"""
Instantaneous :symbols:`radius_of_curvature`.
"""

law = Eq(centripetal_acceleration, speed**2 / radius_of_curvature)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from acceleration and velocity definitions

## Let's assume we are having movement in 2-D space.
## Object position is described with it's radius-vector R - the vector from zero coordinates to the object and with angle '_alpha' between X-axis and this radius-vector.

_time = symbols.time
_alpha = clone_as_function(symbols.angular_distance, [_time])
_cartesian_coordinates = CoordinateSystem()

_curve_radius_horisontal = projector.law.rhs.subs({
    projector.vector_length: radius_of_curvature,
    projector.vector_angle: _alpha(_time)
})
_curve_radius_vertical = projector.law.rhs.subs({
    projector.vector_length: radius_of_curvature,
    projector.vector_angle: pi / 2 - _alpha(_time)
})

## Velocity projections are derivatives of respective coordinates.

#NOTE: replace 'moving_time' first as Derivative can have difficulties when processing both substitutions at once
_velocity_horisontal = velocity_def.definition.rhs.subs(velocity_def.time,
    _time).subs(velocity_def.distance(_time), _curve_radius_horisontal).doit()
_velocity_vertical = velocity_def.definition.rhs.subs(velocity_def.time,
    _time).subs(velocity_def.distance(_time), _curve_radius_vertical).doit()
_velocity_vector = Vector([_velocity_horisontal, _velocity_vertical], _cartesian_coordinates)

## These unit vectors should not necessary be derived. We can choose them at will and prove that
## they are orthogonal to each other and _radial_unit_vector is orthogonal to '_velocity_vector'.
## One can also show that '_tangential_unit_vector' is '_radial_unit_vector' derivative.
_radial_unit_vector = Vector([cos(_alpha(_time)), sin(_alpha(_time))], _cartesian_coordinates)
_tangential_unit_vector = Vector([-sin(_alpha(_time)), cos(_alpha(_time))], _cartesian_coordinates)

## This is Dot product of radial vector and velocity vector. Radial vector is orthogonal to velocity hence vector
## multiplication result should be zero.
assert expr_equals(dot_vectors(_radial_unit_vector, _velocity_vector), 0)
## Radial vector is orthogonal to tangential vector hence tangential vector should be parallel to velocity vector.
assert expr_equals(dot_vectors(_tangential_unit_vector, _radial_unit_vector), 0)

## Use acceleration definition to calculate '_acceleration_vector'
_acceleration_horisontal = acceleration_def.definition.rhs.subs(acceleration_def.time, _time)
_acceleration_horisontal = _acceleration_horisontal.subs(acceleration_def.speed(_time),
    _velocity_horisontal).doit()
_acceleration_vertical = acceleration_def.definition.rhs.subs(acceleration_def.time, _time)
_acceleration_vertical = _acceleration_vertical.subs(acceleration_def.speed(_time),
    _velocity_vertical).doit()
_acceleration_vector = Vector([_acceleration_horisontal, _acceleration_vertical],
    _cartesian_coordinates)

## Prove that '_acceleration_vector' has tangential and radial parts.

_tangential_acceleration_magnitude = radius_of_curvature * Derivative(_alpha(_time), (_time, 2))
_radial_acceleration_magnitude = -radius_of_curvature * Derivative(_alpha(_time), _time)**2

## Use Dot product to find tangential and radial components of acceleration. Confirm they are
## equal to expected value: _tangential_acceleration_magnitude, _radial_acceleration_magnitude
_tangential_acceleration_component = dot_vectors(_acceleration_vector, _tangential_unit_vector)
_radial_acceleration_component = dot_vectors(_acceleration_vector, _radial_unit_vector)
assert expr_equals(_tangential_acceleration_component, _tangential_acceleration_magnitude)
assert expr_equals(_radial_acceleration_component, _radial_acceleration_magnitude)

## Here we've proven that tangential_acceleration + radial_acceleration equals to _acceleration_vector. It means, we've
## changed basis of _acceleration_vector to tangential and radial vectors instead of cartesian coordinates.
## Same result could be achieved by rotating coordinate system by velocity vector angle.

## We are not interested in tangential_acceleration as we are looking for centripetal acceleration which is 'radial_acceleration'
## in our proof.

_angular_velocity_applied = angular_velocity_def.definition.rhs.subs(angular_velocity_def.time,
    _time)
_angular_velocity_applied = _angular_velocity_applied.subs(
    angular_velocity_def.angular_distance(_time), _alpha(_time))
_linear_velocity_applied = linear_velocity_law.law.rhs.subs({
    linear_velocity_law.angular_speed: _angular_velocity_applied,
    linear_velocity_law.radius_of_curvature: radius_of_curvature
})
_law_acceleration = law.rhs.subs(speed, _linear_velocity_applied)

## _radial_acceleration_magnitude has minus sign. It means it is directed towards the center of the curve. The centripetal
## acceleration law is not defined in vector terms so we should only compare acceleration magnitudes (absolute values).
assert expr_equals_abs(_radial_acceleration_magnitude, _law_acceleration)


@validate_input(linear_velocity_=speed, curve_radius_=radius_of_curvature)
@validate_output(centripetal_acceleration)
def calculate_acceleration(linear_velocity_: Quantity, curve_radius_: Quantity) -> Quantity:
    solved = solve(law, centripetal_acceleration, dict=True)[0][centripetal_acceleration]
    result_expr = solved.subs({speed: linear_velocity_, radius_of_curvature: curve_radius_})
    return Quantity(result_expr)
