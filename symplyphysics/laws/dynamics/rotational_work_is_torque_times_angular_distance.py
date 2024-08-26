"""
Rotational work is torque times angular distance
================================================

When a torque accelerates a rigid body in rotation about a fixed axis, the torque does work on
the body. When the torque is constant, the work done on the body is proportional to torque
and the angular displacement during the movement of the body.

..
    TODO Rename file
"""

from sympy import Eq, integrate, solve, pi
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import mechanical_work_from_force_and_move as linear_work_law
from symplyphysics.laws.kinematics import angular_position_is_arc_length_over_radius as angular_position_def
from symplyphysics.laws.dynamics import torque_via_force_and_radial_distance as torque_def
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_law

work = Symbol("work", units.energy)
"""
The work done by the torque.

Symbol:
    :code:`W`
"""

torque = Symbol("torque", units.force * units.length)
r"""
The torque acting on the body.

Symbol:
    :code:`tau`

Latex:
    :math:`\tau`
"""

angular_distance = Symbol("angular_distance", angle_type)
r"""
The angular displacement of the body.

Symbol:
    :code:`delta(theta)`

Latex:
    :math:`\Delta \theta`
"""

law = Eq(work, torque * angular_distance)
r"""
:code:`W = tau * delta(theta)`

Latex:
    .. math::
        W = \tau \Delta \theta
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

force = symbols.dynamics.force
angle = Symbol("angle", angle_type)  # angle between force and radius vectors
radius = Symbol("radius", units.length)

tangent_force = solve(projection_law.law, projection_law.projection)[0].subs({
    projection_law.vector_length: force,
    projection_law.vector_angle: pi / 2 - angle  # complementary angle
})

distance_traveled = solve(
    angular_position_def.law,
    angular_position_def.arc_length,
)[0].subs({
    angular_position_def.angular_position: angular_distance,
    angular_position_def.distance_to_axis: radius,
})

work_derived = linear_work_law.law.rhs.subs({
    linear_work_law.force: tangent_force,
    linear_work_law.distance: distance_traveled,
})

torque_def_sub = torque_def.law.subs({
    torque_def.torque: torque,
    torque_def.force: force,
    torque_def.radial_distance: radius,
    torque_def.angle_between_vectors: angle,
})
work_derived_sub = solve([Eq(work, work_derived), torque_def_sub], (radius, work),
    dict=True)[0][work]

assert expr_equals(work_derived_sub, law.rhs)

# Having derived for the infinitesimal case, let us prove this equality in case of a finite angular
# displacement, with the torque being constant.

angle_start = Symbol("angle_start", angle_type)
angle_end = Symbol("angle_end", angle_type)
# Infinitesimal work = Tau(angle) * dAngle
# SymPy integration does not need dAngle
# And we do not substitute torque for function because it is constant
work_function = work_derived_sub.subs(angular_distance, 1)

integral_work = integrate(work_function, (angle, angle_start, angle_end))

integral_work_solved = solve(
    [Eq(angular_distance, angle_end - angle_start),
    Eq(work, integral_work)], (angle_end, work),
    dict=True)[0][work]

assert expr_equals(integral_work_solved, law.rhs)


@validate_input(torque_=torque, angular_displacement_=angular_distance)
@validate_output(work)
def calculate_work(torque_: Quantity, angular_displacement_: Quantity | float) -> Quantity:
    result = law.rhs.subs({
        torque: torque_,
        angular_distance: angular_displacement_,
    })
    return Quantity(result)
