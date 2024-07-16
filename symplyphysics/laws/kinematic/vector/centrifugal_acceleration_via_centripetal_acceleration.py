"""
Centrifugal acceleration via centripetal acceleration
=====================================================

The vector of centrifugal acceleration has the same magnitude as the vector of centripetal
acceleration but is directed oppositely to it.
"""

from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    scale_vector,
)


def centrifugal_law(centripetal_acceleration_: Vector) -> Vector:
    r"""Centrifugal acceleration via centripetal acceleration.

    Law:
        :math:`\vec a_\text{cf} = - \vec a_\text{cf}`

    @param centripetal_acceleration_ The acceleration experienced by a rotating body in an inertial frame.

    @return The acceleration experienced by a body in a non-inertial, rotating frame.
    """
    return scale_vector(-1, centripetal_acceleration_)


def centripetal_law(centrifugal_acceleration_: Vector) -> Vector:
    r"""Centripetal acceleration via centrifugal acceleration.

    Law:
        :math:`\vec a_\text{cp} = - \vec a_\text{cf}`
    
    @param centrifugal_acceleration_ The acceleration experienced by a body in a non-inertial, rotating frame.

    @return The acceleration experienced by a rotating body in an inertial frame.
    """

    return scale_vector(-1, centrifugal_acceleration_)


@validate_input(centripetal_acceleration_=units.acceleration)
@validate_output(units.acceleration)
def calculate_centrifugal_acceleration(
    centripetal_acceleration_: QuantityVector,
) -> QuantityVector:
    result_vector = centrifugal_law(centripetal_acceleration_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector)
