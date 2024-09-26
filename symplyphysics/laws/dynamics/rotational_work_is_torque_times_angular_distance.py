"""
Rotational work is torque times angular distance
================================================

When a torque accelerates a rigid body in rotation about a fixed axis, the torque does work on
the body. When the torque is constant, the work done on the body is proportional to torque
and the angular displacement during the movement of the body.
"""

from sympy import Eq, integrate, solve, pi
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import mechanical_work_from_force_and_distance as linear_work_law
from symplyphysics.laws.kinematics import angular_position_is_arc_length_over_radius as angular_position_def
from symplyphysics.laws.dynamics import torque_via_force_and_radial_distance as torque_def
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_law

work = symbols.work
"""
The :symbols:`work` done by the torque.
"""

torque = symbols.torque
r"""
The :symbols:`torque` acting on the body.
"""

angular_distance = symbols.angular_distance
r"""
The :symbols:`angular_distance` of the body.
"""

law = Eq(work, torque * angular_distance)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from non-rotational counterpart.
## We assume a particle moving along a curved trajectory due to a force F applied to it.
## Since only the component of the force vector which accelerates the particle does work on it,
## we are only interested in the tangent component of it.

# Conditions:
## The displacement of the particle should be small enough, so that the force, radius vector,
## and angle between the two stay the same at all points of the path in question.

# Reference frame:
## the x axis points in the direction of the radius vector
## the y axis is tangent to the part of the path in question

_force = symbols.force
_angle = symbols.angle  # angle between force and radius vectors
_radius = symbols.distance_to_axis

_tangent_force = solve(projection_law.law, projection_law.projection)[0].subs({
    projection_law.vector_length: _force,
    projection_law.vector_angle: pi / 2 - _angle  # complementary angle
})

_distance_traveled = solve(
    angular_position_def.law,
    angular_position_def.arc_length,
)[0].subs({
    angular_position_def.angular_position: angular_distance,
    angular_position_def.distance_to_axis: _radius,
})

_work_derived = linear_work_law.law.rhs.subs({
    linear_work_law.force: _tangent_force,
    linear_work_law.distance: _distance_traveled,
})

_torque_def_sub = torque_def.law.subs({
    torque_def.torque: torque,
    torque_def.force: _force,
    torque_def.radial_distance: _radius,
    torque_def.angle_between_vectors: _angle,
})
_work_derived_sub = solve([Eq(work, _work_derived), _torque_def_sub], (_radius, work),
    dict=True)[0][work]

assert expr_equals(_work_derived_sub, law.rhs)

# Having derived for the infinitesimal case, let us prove this equality in case of a finite angular
# displacement, with the torque being constant.

_angle_start = clone_as_symbol(symbols.angle, display_symbol="phi_0")
_angle_end = clone_as_symbol(symbols.angle, display_symbol="phi_1")
# Infinitesimal work = Tau(angle) * dAngle
# SymPy integration does not need dAngle
# And we do not substitute torque for function because it is constant
_work_function = _work_derived_sub.subs(angular_distance, 1)

_integral_work = integrate(_work_function, (_angle, _angle_start, _angle_end))

_integral_work_solved = solve(
    [Eq(angular_distance, _angle_end - _angle_start),
    Eq(work, _integral_work)], (_angle_end, work),
    dict=True)[0][work]

assert expr_equals(_integral_work_solved, law.rhs)


@validate_input(torque_=torque, angular_displacement_=angular_distance)
@validate_output(work)
def calculate_work(torque_: Quantity, angular_displacement_: Quantity | float) -> Quantity:
    result = law.rhs.subs({
        torque: torque_,
        angular_distance: angular_displacement_,
    })
    return Quantity(result)
