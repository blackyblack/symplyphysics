r"""
Centrifugal acceleration via centripetal acceleration
=====================================================

Centrifugal acceleration has the same magnitude as centripetal acceleration but is directed oppositely to it.

Law
===
    :math:`\vec a_\text{cf} = -\vec a_\text{cp}`
"""

from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    scale_vector,
)

centripetal_acceleration_: Vector
r"""
The acceleration experienced by a rotating body in an inertial frame

Symbol:
    :math:`\vec a_\text{cp}`
"""

centrifugal_acceleration_: Vector
r"""
The acceleration experienced by a body in a non-inertial, rotating frame

Symbol:
    :math:`\vec a_\text{cf}`
"""


def centrifugal_law(centripetal_acceleration_: Vector) -> Vector:
    return scale_vector(-1, centripetal_acceleration_)


def centripetal_law(centrifugal_acceleration_: Vector) -> Vector:
    return scale_vector(-1, centrifugal_acceleration_)


@validate_input(centripetal_acceleration_=units.acceleration)
@validate_output(units.acceleration)
def calculate_centrifugal_acceleration(
    centripetal_acceleration_: QuantityVector,
) -> QuantityVector:
    result_vector = centrifugal_law(centripetal_acceleration_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector)
