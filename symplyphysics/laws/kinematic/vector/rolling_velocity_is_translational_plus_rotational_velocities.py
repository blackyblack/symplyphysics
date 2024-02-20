from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    add_cartesian_vectors,
)

# Description
## Suppose a wheel is rolling smoothly along a surface, i.e. without slipping or bouncing on the
## surface. Such motion can be viewed as a combination of translational motion of the wheel's center
## of mass, and rotational motion of its points about the rotational axis.

# Law: v = v_translation + v_rotation
## v - velocity of a point on the wheel relative to the inertial frame
## v_translation - velocity of the wheel's center of mass relative to the inertial frame
## v_rotation - see [rotational velocity formula](./linear_velocity_is_angular_velocity_cross_radius.py)
##              relative to wheel


def rolling_velocity_law(translational_velocity_: Vector, rotational_velocity_: Vector) -> Vector:
    return add_cartesian_vectors(translational_velocity_, rotational_velocity_)


@validate_input(translational_velocity_=units.velocity, rotational_velocity_=units.velocity)
@validate_output(units.velocity,)
def calculate_rolling_velocity(
    translational_velocity_: QuantityVector, rotational_velocity_: QuantityVector,
) -> QuantityVector:
    result = rolling_velocity_law(
        translational_velocity_.to_base_vector(),
        rotational_velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(result)
