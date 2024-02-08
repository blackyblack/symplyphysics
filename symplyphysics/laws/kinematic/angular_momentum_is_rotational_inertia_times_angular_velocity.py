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
## For a rigid body rotating about a fixed axis, the component of its angular momentum parallel
## to the rotational axis is found as the product of the body's rotational inertia and the magnitude
## of its angular velocity.

# Law: L = I * w
## L - component of angular momentum parallel to rotational axis
## I - rotational inertia
## w - angular velocity

angular_momentum = Symbol("angular_momentum", units.length * units.momentum)
rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
angular_velocity = Symbol("angular_velocity", angle_type / units.time)

law = Eq(angular_momentum, rotational_inertia * angular_velocity)


def print_law() -> str:
    return print_expression(law)


@validate_input(rotational_inertia_=rotational_inertia, angular_velocity_=angular_velocity)
@validate_output(angular_momentum)
def calculate_angular_momentum(rotational_inertia_: Quantity, angular_velocity_: Quantity) -> Quantity:
    result = law.rhs.subs({
        rotational_inertia: rotational_inertia_,
        angular_velocity: angular_velocity_,
    })
    return Quantity(result)
