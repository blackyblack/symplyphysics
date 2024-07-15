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
## Suppose two reference frames, one of which is fixed (S) and the other one is moving (S'). When
## the body is moving within a rotating coordinate system, its path deflects due to the appearance
## of the Coriolis acceleration on it. The object does not actually deviate from its path per se
## but it appears to do so because of the motion of the coordinate system.

## Suppose a reference frame S' is fixed to a rotating body A (e.g. Earth), so that frame S' rotates w.r.t.
## another static reference frame S. The Coriolis acceleration is the acceleration another body B has when
## moving within rotating reference frame S', so it is essentially zero for objects at rest in S'.

# Law: a_cor = 2 * cross(v_rel, w)
## a_cor - vector of Coriolis acceleration of body B in S'
## v_rel - vector of velocity of body B in moving frame S'
## w - pseudovector of angular velocity of rotation of moving frame S'


def coriolis_acceleration_law(
    velocity_: Vector,
    angular_velocity_: Vector,
) -> Vector:
    return scale_vector(2, cross_cartesian_vectors(velocity_, angular_velocity_))


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
        velocity_.to_base_vector(),
        angular_velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(result)
