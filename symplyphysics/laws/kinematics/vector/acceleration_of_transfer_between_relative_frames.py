r"""
Acceleration of transfer between relative frames
================================================

Imagine two reference frames, one of which is fixed (:math:`S`) and the other one is moving (:math:`S'`). The motion of
a body stationary in moving frame :math:`S'` due to the motion of the frame itself is called transfer motion.
The acceleration related to such motion is called transfer acceleration. It is composed of the acceleration
of the moving frame relative to the fixed frame, centripetal acceleration and the acceleration due to uneven
rotation of the moving frame. The transfer acceleration only depends on the motion of frame :math:`S'` relative to
stationary frame :math:`S`, so its physical meaning would be that it is the acceleration in :math:`S` of a point stationary
in :math:`S'`.
"""

from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    add_cartesian_vectors,
    subtract_cartesian_vectors,
)

def transfer_acceleration_law(
    moving_frame_acceleration_: Vector,
    centripetal_acceleration_: Vector,
    rotation_acceleration_: Vector,
) -> Vector:
    r"""
    Transfer acceleration as a sum of accelerations.

    Law:
        :code:`a_tr = a_0 + a_cp + a_rot`

    Latex:
        .. math::
            {\vec a}_\text{tr} = {\vec a}_0 + {\vec a}_\text{cp} + {\vec a}_\text{rot}

    :param moving\_frame\_acceleration\_: acceleration of :math:`S'` relative to :math:`S`

        Symbol: :code:`a_0`

        Latex: :math:`{\vec a}_0`

        Dimension: *acceleration*

    :param centripetal\_acceleration\_: centripetal acceleration of body in :math:`S'`

        Symbol: :code:`a_cp`

        Latex: :math:`{\vec a}_\text{cp}`

        Dimension: *acceleration*

    :param rotation\_acceleration\_: acceleration caused by non-uniform rotation of :math:`S'`

        Symbol: :code:`a_rot`

        Latex: :math:`{\vec a}_\text{rot}`

        Dimension: *acceleration*

    :return: transfer acceleration of body

        Symbol: :code:`a_tr`

        Latex: :math:`{\vec a}_\text{tr}`

        Dimension: *acceleration*
    """

    return add_cartesian_vectors(
        moving_frame_acceleration_,
        add_cartesian_vectors(
        centripetal_acceleration_,
        rotation_acceleration_,
        ))


# a_0 = a_tr - (a_centripetal + a_non_uniform_rotation)
def moving_frame_acceleration_law(
    transfer_acceleration_: Vector,
    centripetal_acceleration_: Vector,
    rotation_acceleration_: Vector,
) -> Vector:
    r"""
    Acceleration of :math:`S'` relative to :math:`S`.

    Law:
        :code:`a_0 = a_tr - (a_cp + a_rot)`

    Latex:
        .. math::
            {\vec a}_0 = {\vec a}_\text{tr} - ({\vec a}_\text{cp} + {\vec a}_\text{rot})

    :param transfer\_acceleration\_: transfer acceleration of body

        Symbol: :code:`a_tr`

        Latex: :math:`{\vec a}_\text{tr}`

        Dimension: *acceleration*

    :param centripetal\_acceleration\_: centripetal acceleration of body in :math:`S'`

        Symbol: :code:`a_cp`

        Latex: :math:`{\vec a}_\text{cp}`

        Dimension: *acceleration*

    :param rotation\_acceleration\_: acceleration caused by non-uniform rotation of :math:`S'`

        Symbol: :code:`a_rot`

        Latex: :math:`{\vec a}_\text{rot}`

        Dimension: *acceleration*
    
    :return: acceleration of :math:`S'` relative to :math:`S`

        Symbol: :code:`a_0`

        Latex: :math:`{\vec a}_0`

        Dimension: *acceleration*
    """

    return subtract_cartesian_vectors(
        transfer_acceleration_,
        add_cartesian_vectors(
        centripetal_acceleration_,
        rotation_acceleration_,
        ))


# a_centripetal = a_tr - (a_0 + a_non_uniform_rotation)
def centripetal_acceleration_law(
    transfer_acceleration_: Vector,
    moving_frame_acceleration_: Vector,
    rotation_acceleration_: Vector,
) -> Vector:
    r"""
    Centripetal acceleration in :math:`S'`.

    Law:
        :code:`a_cp = a_tr - (a_0 + a_rot)`

    Latex:
        .. math::
            {\vec a}_\text{cp} = {\vec a}_\text{tr} - ({\vec a}_0 + {\vec a}_\text{rot})

    :param transfer\_acceleration\_: transfer acceleration of body

        Symbol: :code:`a_tr`

        Latex: :math:`{\vec a}_\text{tr}`

        Dimension: *acceleration*

    :param moving\_frame\_acceleration\_: acceleration of :math:`S'` relative to :math:`S`

        Symbol: :code:`a_0`

        Latex: :math:`{\vec a}_0`

        Dimension: *acceleration*

    :param rotation\_acceleration\_: acceleration caused by non-uniform rotation of :math:`S'`

        Symbol: :code:`a_rot`

        Latex: :math:`{\vec a}_\text{rot}`

        Dimension: *acceleration*

    :return: centripetal acceleration of body in :math:`S'`

        Symbol: :code:`a_cp`

        Latex: :math:`{\vec a}_\text{cp}`

        Dimension: *acceleration*
    """

    return subtract_cartesian_vectors(
        transfer_acceleration_,
        add_cartesian_vectors(
        moving_frame_acceleration_,
        rotation_acceleration_,
        ))


# a_non_uniform_rotation = a_tr - (a_0 + a_centripetal)
def rotation_acceleration_law(
    transfer_acceleration_: Vector,
    moving_frame_acceleration_: Vector,
    centripetal_acceleration_: Vector,
) -> Vector:
    r"""
    Acceleration due to non-uniform rotation of :math:`S'`.

    Law:
        :code:`a_rot = a_tr - (a_0 + a_cp)`

    Latex:
        .. math::
            {\vec a}_\text{rot} = {\vec a}_\text{tr} - ({\vec a}_0 + {\vec a}_\text{cp})

    :param transfer\_acceleration\_: transfer acceleration of body

        Symbol: :code:`a_tr`

        Latex: :math:`{\vec a}_\text{tr}`

        Dimension: *acceleration*

    :param moving\_frame\_acceleration\_: acceleration of :math:`S'` relative to :math:`S`

        Symbol: :code:`a_0`

        Latex: :math:`{\vec a}_0`

        Dimension: *acceleration*

    :param centripetal\_acceleration\_: centripetal acceleration of body in :math:`S'`

        Symbol: :code:`a_cp`

        Latex: :math:`{\vec a}_\text{cp}`

        Dimension: *acceleration*

    :return: acceleration caused by non-uniform rotation of :math:`S'`

        Symbol: :code:`a_rot`

        Latex: :math:`{\vec a}_\text{rot}`

        Dimension: *acceleration*
    """

    return subtract_cartesian_vectors(
        transfer_acceleration_,
        add_cartesian_vectors(
        moving_frame_acceleration_,
        centripetal_acceleration_,
        ))


@validate_input(
    moving_frame_acceleration_=units.acceleration,
    centripetal_acceleration_=units.acceleration,
    rotation_acceleration_=units.acceleration,
)
@validate_output(units.acceleration)
def calculate_transfer_acceleration(
    moving_frame_acceleration_: QuantityVector,
    centripetal_acceleration_: QuantityVector,
    rotation_acceleration_: QuantityVector,
) -> QuantityVector:
    vector = transfer_acceleration_law(
        moving_frame_acceleration_.to_base_vector(),
        centripetal_acceleration_.to_base_vector(),
        rotation_acceleration_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(vector)
