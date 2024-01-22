from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)

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
