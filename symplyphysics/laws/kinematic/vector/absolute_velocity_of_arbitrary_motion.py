from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    add_cartesian_vectors,
    scale_vector,
)

# Description
## Suppose two reference frames, one of which is fixed (S) and the other one is moving arbitrarily (S'). The movement of the
## body relative to fixed frame S is called absolute movement. The movement of the body relative to moving frame S'
## is called relative movement. And the movement of the body due to the movement of reference frame S' is called transfer
## movement. The absolute velocity is the sum of relative and transfer velocities.

# Law: v_abs = v_rel + v_tr
## v_abs - vector of absolute velocity relative to fixed frame S
## v_rel - vector of velocity relative to moving frame S'
## v_tr - vector of velocity of transfer between frames S and S'

# Notes
## - Moving frame S' can perform both translational and rotational motion.


def absolute_velocity_law(
    relative_velocity_: Vector,
    transfer_velocity_: Vector,
) -> Vector:
    return add_cartesian_vectors(
        relative_velocity_,
        transfer_velocity_,
    )


def relative_velocity_law(
    absolute_velocity_: Vector,
    transfer_velocity_: Vector,
) -> Vector:
    return add_cartesian_vectors(absolute_velocity_, scale_vector(-1, transfer_velocity_))


def transfer_velocity_law(
    absolute_velocity_: Vector,
    relative_velocity_: Vector,
) -> Vector:
    return add_cartesian_vectors(absolute_velocity_, scale_vector(-1, relative_velocity_))


@validate_input(
    relative_velocity_=units.velocity,
    transfer_velocity_=units.velocity,
)
@validate_output(units.velocity)
def calculate_absolute_velocity(
    relative_velocity_: QuantityVector,
    transfer_velocity_: QuantityVector,
) -> QuantityVector:
    result = absolute_velocity_law(
        relative_velocity_.to_base_vector(),
        transfer_velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(result)
