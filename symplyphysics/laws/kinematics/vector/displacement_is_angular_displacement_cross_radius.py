"""
Linear displacement is angular displacement cross radius
=================================================

Assuming a body rotating around a fixed axis, the vector of its linear displacement can be expressed
as the cross product of the pseudovector of angular displacement and the radius vector of rotation.

**Conditions:**

#. The axis is fixed.
#. Angular displacement pseudovector and radius vector must be orthogonal to one another.
"""

from pytest import approx
from symplyphysics import (
    Quantity,
    dot_vectors,
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    cross_cartesian_vectors,
)


def displacement_law(angular_displacement_: Vector, rotation_radius_: Vector) -> Vector:
    r"""
    Displacement vector.

    Law:
        :code:`s = cross(theta, r)`

    Latex:
        .. math::
            \vec s = \vec \theta \times \vec r

    :param angular_displacement\_: pseudovector of angular displacement parallel to axis of rotation

        Symbol: :code:`theta`
        
        Latex: :math:`\vec \theta`

        Dimension: *angle*

    :param rotation_radius\_: radius vector pointing away from the rotational axis and perpendicular to it

        Symbol: :code:`r`

        Latex: :math:`\vec r`

        Dimension: *length*

    :return: vector of linear displacement

        Symbol: :code:`s`

        Latex: :math:`\vec s`

        Dimension: *length*
    """

    return cross_cartesian_vectors(angular_displacement_, rotation_radius_)


@validate_input(angular_displacement_=angle_type, rotation_radius_=units.length)
@validate_output(units.length)
def calculate_linear_displacement(angular_displacement_: QuantityVector,
    rotation_radius_: QuantityVector) -> QuantityVector:
    angular_displacement_vector = angular_displacement_.to_base_vector()
    rotation_radius_vector = rotation_radius_.to_base_vector()
    dot_vectors_result = Quantity(dot_vectors(angular_displacement_vector, rotation_radius_vector))
    if dot_vectors_result.scale_factor != approx(0.0, rel=1e-3):
        raise ValueError(
            "Angular displacement pseudovector and rotation radius vector should be perpendicular to each other"
        )
    result_vector = displacement_law(angular_displacement_vector, rotation_radius_vector)
    return QuantityVector.from_base_vector(result_vector)
