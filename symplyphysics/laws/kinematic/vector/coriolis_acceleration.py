from symplyphysics import (
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    scale_vector,
    cross_cartesian_vectors,
)

# Description
## Suppose two reference frames, one of which is fixed (S) and the other one is moving (S').
## TODO

# Law: a_cor = 2 * cross(w, v_rel)
## a_cor - vector of Coriolis acceleration
## w - pseudovector of angular velocity
## v_rel - vector of velocity relative to moving frame S'


def coriolis_acceleration_law(
    angular_velocity_: Vector,
    velocity_: Vector,
) -> Vector:
    return scale_vector(2, cross_cartesian_vectors(angular_velocity_, velocity_))


@validate_input(
    angular_velocity_=angle_type / units.time,
    velocity_=units.velocity,
)
@validate_output(units.acceleration)
def calculate_coriolis_acceleration(
    angular_velocity_: QuantityVector,
    velocity_: QuantityVector,
) -> QuantityVector:
    result = coriolis_acceleration_law(
        angular_velocity_.to_base_vector(),
        velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(result)
