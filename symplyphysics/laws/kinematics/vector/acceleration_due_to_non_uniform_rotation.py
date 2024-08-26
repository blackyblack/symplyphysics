r"""
Acceleration due to non-uniform rotation
========================================

Imagine two reference frames, one of which is fixed (:math:`S`) and the other is moving (:math:`S'`). When :math:`S'` rotates
around :math:`S` in a non-uniform way, the acceleration of some body :math:`B` in :math:`S` has a component corresponding to
that non-uniform rotation of :math:`S'`. It is part of the transfer acceleration of body :math:`B` in :math:`S`.

**Notation:**

#. :math:`\vec a \times \vec b` (:code:`cross(a, b)`) is vector product of :math:`\vec a` and :math:`\vec b`.
"""

from sympy import Expr, abc
from symplyphysics import (
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    Quantity,
    QuantityVector,
    cross_cartesian_vectors,
    scale_vector,
)
from symplyphysics.core.vectors.arithmetics import diff_cartesian_vector

def non_uniform_rotation_acceleration_law(
    angular_velocity_: Vector,
    time_: Expr,
    radius_vector_: Vector,
) -> Vector:
    r"""
    Acceleration due to non-uniform rotation.

    Law:
        :code:`a_rot = cross(Derivative(w(t), t), r)`

    Latex:
        .. math::
            {\vec a}_\text{rot} = \frac{d \vec \omega}{d t} \times \vec r

    :param angular_velocity\_: angular velocity as a function of time

        Symbol: :code:`w(t)`

        Latex: :math:`\vec \omega(t)`

        Dimension: *angle* / *time*

    :param time\_: time

        Symbol: :code:`t`

        Dimension: *time*

    :param radius_vector\_: radius vector, or position vector, of body

        Symbol: :code:`r`

        Latex: :math:`\vec r`

        Dimension: *length*

    :return: acceleration due to non-uniform rotation

        Symbol: :code:`a_rot`

        Latex: :math:`{\vec a}_\text{rot}`

        Dimension: *acceleration*
    """

    return cross_cartesian_vectors(
        diff_cartesian_vector(angular_velocity_, time_),
        radius_vector_,
    )


@validate_input(
    angular_velocity_change_=angle_type / units.time,
    time_change_=units.time,
    position_vector=units.length,
)
@validate_output(units.acceleration)
def calculate_non_uniform_rotation_acceleration(
    angular_velocity_change_: QuantityVector,
    time_change_: Quantity,
    radius_vector_: QuantityVector,
) -> QuantityVector:
    time_ = abc.t

    angular_velocity_ = scale_vector(
        time_ / time_change_,
        angular_velocity_change_.to_base_vector(),
    )

    result_vector = non_uniform_rotation_acceleration_law(
        angular_velocity_,
        time_,
        radius_vector_.to_base_vector(),
    )

    return QuantityVector.from_base_vector(result_vector)
