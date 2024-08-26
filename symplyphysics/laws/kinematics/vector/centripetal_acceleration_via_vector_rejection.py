r"""
Centripetal acceleration via vector rejection
=============================================

*Centripetal acceleration* is the acceleration of a body in a rotating coordinate system
which is directed towards the axis of rotation.

Also see :doc:`laws.kinematics.vector.centripetal_acceleration_via_cross_product`.

**Notation:**

#. :math:`|\vec a|` (:code:`norm(a)`) is the Euclidean norm of :math:`\vec a`.
#. :math:`\text{oproj}_{\vec b} \vec a` (:code:`reject(a, b)`) is the rejection of
   :math:`\vec a` from :math:`\vec b`, i.e. the component of :math:`\vec a` orthogonal
   to :math:`\vec b`.
"""

from symplyphysics import (
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    vector_magnitude,
    scale_vector,
)
from symplyphysics.core.vectors.arithmetics import reject_cartesian_vector


def centripetal_acceleration_law(
    angular_velocity_: Vector,
    radius_vector_: Vector,
) -> Vector:
    r"""
    Centripetal acceleration via angular velocity and radius vector.

    Law:
        :code:`a_c = -1 * norm(w)^2 * reject(r, w)`

    Latex:
        .. math::
            {\vec a}_\text{c} = - |\vec \omega|^2 \text{oproj}_{\vec \omega} \vec r

    :param angular_velocity\_: pseudovector of angular velocity

        Symbol: :code:`w`

        Latex: :math:`\vec \omega`

        Dimension: *angle* / *time*

    :param radius_vector\_: radius vector, or position vector

        Symbol: :code:`r`

        Latex: :math:`\vec r`

        Dimension: *length*

    :return: vector of centripetal acceleration

        Symbol: :code:`a_c`

        Latex: :math:`{\vec a}_\text{c}`

        Dimension: *acceleration*
    """

    return scale_vector(
        -1 * vector_magnitude(angular_velocity_)**2,
        reject_cartesian_vector(radius_vector_, angular_velocity_),
    )


@validate_input(
    angular_velocity_=angle_type / units.time,
    radius_vector_=units.length,
)
@validate_output(units.acceleration)
def calculate_centripetal_acceleration(
    angular_velocity_: QuantityVector,
    radius_vector_: QuantityVector,
) -> QuantityVector:
    vector = centripetal_acceleration_law(
        angular_velocity_.to_base_vector(),
        radius_vector_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(vector)
