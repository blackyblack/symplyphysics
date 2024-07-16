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
    r"""
    Centrifugal acceleration via centripetal acceleration.

    Law:
        a_centrifugal = -a_centripetal

    Latex:
        :math:`\vec a_\text{cf} = - \vec a_\text{cp}`

    :param centripetal_acceleration\_: The :attr:`~symplyphysics.symbols.kinematic.acceleration` experienced by a rotating body in an inertial frame.

        Dimension: *acceleration*
    
    :return: The acceleration experienced by a body in a non-inertial, rotating frame.

        Dimension: *acceleration*
    """
    return scale_vector(-1, centripetal_acceleration_)


def centripetal_law(centrifugal_acceleration_: Vector) -> Vector:
    r"""
    Centripetal acceleration via centrifugal acceleration.

    Law:
        a_centripetal = -a_centrifugal
    
    Latex:
        :math:`\vec a_\text{cp} = - \vec a_\text{cf}`
    
    :param centrifugal_acceleration\_: The :attr:`~symplyphysics.symbols.kinematic.acceleration` experienced by a body in a non-inertial, rotating frame.
    
        Dimension: *acceleration*

    :return: The :attr:`~symplyphysics.symbols.kinematic.acceleration` experienced by a rotating body in an inertial frame.
    
        Dimension: *acceleration*
    """

    return scale_vector(-1, centrifugal_acceleration_)


@validate_input(centripetal_acceleration_=units.acceleration)
@validate_output(units.acceleration)
def calculate_centrifugal_acceleration(
    centripetal_acceleration_: QuantityVector,
) -> QuantityVector:
    result_vector = centrifugal_law(centripetal_acceleration_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector)
