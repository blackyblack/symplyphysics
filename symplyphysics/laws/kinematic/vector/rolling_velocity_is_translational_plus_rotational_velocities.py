from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    add_cartesian_vectors,
)

# Description
## Suppose an object is rolling smoothly along a surface, i.e. without slipping or bouncing on the
## surface. Such motion can be viewed as a combination of translational motion of the object's center
## of mass, and rotational motion of its points about the rotational axis.

# Law: v = v_translation + v_rotation
## v - velocity vector of a point on the object.
## v_translation - velocity component responsible for translational movement of object's center
## v_rotation - velocity component responsible for rotational movement of the point on object,
##              see [rotational velocity formula](./linear_velocity_is_angular_velocity_cross_radius.py)

# Notes
## Suppose the object in question is a wheel of constant radius, then
## - v_translation is a vector related to the movement of the wheel's center
##   its magnitude is given by [linear velocity formula](../linear_velocity_from_angular_velocity_and_radius)
##   where radius is the total radius of the wheel
## - v_rotation is calculated via [rotational velocity formula](./linear_velocity_is_angular_velocity_cross_radius.py)
##   where rotation_radius vector is the vector from the wheel's center to the point in question


def rolling_velocity_law(translational_velocity_: Vector, rotational_velocity_: Vector) -> Vector:
    return add_cartesian_vectors(translational_velocity_, rotational_velocity_)


@validate_input(translational_velocity_=units.velocity, rotational_velocity_=units.velocity)
@validate_output(units.velocity,)
def calculate_rolling_velocity(
    translational_velocity_: QuantityVector, rotational_velocity_: QuantityVector,
) -> QuantityVector:
    result = rolling_velocity_law(translational_velocity_, rotational_velocity_)
    return QuantityVector(result.components, translational_velocity_.coordinate_system)
