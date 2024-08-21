"""
Centripetal acceleration via cross product
==========================================

*Centripetal acceleration* is the acceleration of a body in a rotating coordinate system
which is directed towards the axis of rotation.
"""

from sympy import symbols
from symplyphysics import (
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    cross_cartesian_vectors,
    subtract_cartesian_vectors,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics.vector import centripetal_acceleration_via_vector_rejection as rejection_law


def centripetal_acceleration_law(
    angular_velocity_: Vector,
    position_vector_: Vector,
) -> Vector:
    r"""
    Centripetal acceleration via angular velocity and radius vector.

    Law:
        :code:`a_c = cross(w, cross(w, r))`

    Latex:
        .. math::
            {\vec a}_\text{c} = \vec \omega \times (\vec \omega \times \vec r)

    :param angular_velocity\_: pseudovector of angular velocity

        Symbol: :code:`w`

        Latex: :math:`\vec \omega`

        Dimension: *angle* / *time*

    :param position_vector\_: radius vector, or position vector

        Symbol: :code:`r`

        Latex: :math:`\vec r`

        Dimension: *length*

    :return: vector of centripetal acceleration

        Symbol: :code:`a_c`

        Latex: :math:`{\vec a}_\text{c}`

        Dimension: *acceleration*
    """

    return cross_cartesian_vectors(
        angular_velocity_,
        cross_cartesian_vectors(angular_velocity_, position_vector_),
    )


# Prove that these two forms of the law are equivalent.

_angular_velocity = Vector(symbols("angular_velocity_x:z"))
_position_vector = Vector(symbols("position_vector_x:z"))
_cross_product_result = centripetal_acceleration_law(_angular_velocity, _position_vector)
_rejection_result = rejection_law.centripetal_acceleration_law(_angular_velocity, _position_vector)
_difference = subtract_cartesian_vectors(_cross_product_result, _rejection_result).simplify()
for _component in _difference.components:
    expr_equals(_component, 0)


@validate_input(
    angular_velocity_=angle_type / units.time,
    position_vector_=units.length,
)
@validate_output(units.acceleration)
def calculate_centripetal_acceleration(
    angular_velocity_: QuantityVector,
    position_vector_: QuantityVector,
) -> QuantityVector:
    vector = centripetal_acceleration_law(
        angular_velocity_.to_base_vector(),
        position_vector_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(vector)
