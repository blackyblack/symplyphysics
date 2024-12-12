r"""
Coriolis acceleration
=====================

Imagine two reference frames, one of which is fixed (:math:`S`) and the other one is moving (:math:`S'`). When
the body is moving within a rotating coordinate system, its path deflects due to the appearance
of the Coriolis acceleration on it. The object does not actually deviate from its path per se
but it appears to do so because of the motion of the coordinate system.

Suppose a reference frame :math:`S'` is fixed to a rotating body :math:`A` (e.g. Earth), so that frame :math:`S'` rotates w.r.t.
another static reference frame :math:`S`. The Coriolis acceleration is the acceleration another body :math:`B` has when
moving within rotating reference frame :math:`S'`, so it is essentially zero for objects at rest in :math:`S'`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Coriolis_force#Formula>`__.
"""

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


def coriolis_acceleration_law(
    velocity_: Vector,
    angular_velocity_: Vector,
) -> Vector:
    r"""
    Coriolis acceleration via relative velocity of body and angular velocity of rotation of :math:`S'`.

    Law:
        :code:`a_cor = 2 * cross(v_rel, w)`

    Latex:
        .. math::
            {\vec a}_\text{cor} = 2 ({\vec v}_\text{rel} \times \vec \omega)

    :param velocity\_: velocity of body relative to :math:`S'`

        Symbol: :code:`v_rel`

        Latex: :math:`{\vec v}_\text{rel}`

        Dimension: *velocity*

    :param angular\_velocity\_: angular velocity of rotation of :math:`S'` relative to :math:`S`

        Symbol: :code:`w`

        Latex: :math:`\vec \omega`

        Dimension: *angle* / *time*

    :return: Coriolis acceleration of body in :math:`S'`

        Symbol: :code:`a_cor`

        Latex: :math:`{\vec a}_\text{cor}`

        Dimension: *acceleration*
    """

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
