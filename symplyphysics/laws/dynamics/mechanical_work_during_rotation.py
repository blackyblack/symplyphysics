from sympy import Eq, sin, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import mechanical_work_from_force_and_move as linear_work_law
from symplyphysics.laws.kinematic import angular_position_is_arc_length_over_radius as angular_position_def
from symplyphysics.laws.dynamics import torque_due_to_twisting_force as torque_def

# Description
## When a torque accelerates a rigid body in rotation about a fixed axis, the torque does work on
## the body. When the torque is constant, the work done on the body is proportional to torque
## and the angular displacement during the movement of the body.

# Law: W = tau * delta_theta
## W - work done
## tau - torque accelerating the body
## delta_theta - angular displacement of the body

work = Symbol("work", units.energy)
torque = Symbol("torque", units.force * units.length)
angular_displacement = Symbol("angular_displacement", angle_type)

law = Eq(work, torque * angular_displacement)


# Derive law from non-rotational counterpart.
## We assume a particle moving along a curved trajectory due to a force F applied to it.
## The displacement of the particle should be small enough so that the position vector
## is approximately the same at all points of the trajectory. In that case, the force vector
## is assumed to form an angle phi with the position vector of the particle. Since only the
## component of the force vector that is tangent to the path is the one accelerating the
## particle, it does work on the particle.

force = Symbol("force", units.force)
angle = Symbol("angle", angle_type)
radius = Symbol("radius", units.length)

tangent_force = force * sin(angle)

distance_traveled = solve(
    angular_position_def.law,
    angular_position_def.arc_length,
)[0].subs({
    angular_position_def.angular_position: angular_displacement,
    angular_position_def.path_radius: radius,
})

work_derived = linear_work_law.law.rhs.subs({
    linear_work_law.force: tangent_force,
    linear_work_law.distance: distance_traveled,
})

torque_def_sub = torque_def.law.subs({
    torque_def.torque: torque,
    torque_def.force: force,
    torque_def.distance_to_axis: radius,
    torque_def.angle: angle,
})
work_derived_sub = solve(
    [Eq(work, work_derived), torque_def_sub],
    (radius, work),
    dict=True
)[0][work]

work_from_law = law.rhs

assert expr_equals(work_derived_sub, work_from_law)


def print_law() -> str:
    return print_expression(law)


@validate_input(torque_=torque, angular_displacement_=angular_displacement)
@validate_output(work)
def calculate_work(torque_: Quantity, angular_displacement_: Quantity | float) -> Quantity:
    result = law.rhs.subs({
        torque: torque_,
        angular_displacement: angular_displacement_,
    })
    return Quantity(result)
