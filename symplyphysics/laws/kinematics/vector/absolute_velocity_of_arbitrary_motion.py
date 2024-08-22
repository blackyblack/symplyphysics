r"""
Absolute velocity of arbitrary motion
=====================================

Imagine two reference frames, one of which is fixed (:math:`S`) and the other one is moving arbitrarily (:math:`S'`). The motion of the
body relative to fixed frame :math:`S` is called *absolute motion*. The motion of the body relative to moving frame :math:`S'`
is called *relative motion*. The motion of the body due to the motion of reference frame :math:`S'` is called *transfer
motion*. Absolute velocity is the sum of relative and transfer velocities.

**Notes:**

#. Moving frame :math:`S'` can perform both translational and rotational motion.
"""

from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    add_cartesian_vectors,
    scale_vector,
)

# Law: v_abs = v_rel + v_tr
## v_abs - vector of absolute velocity relative to fixed frame S
## v_rel - vector of velocity relative to moving frame S'
## v_tr - vector of velocity of transfer between frames S and S'

def absolute_velocity_law(
    relative_velocity_: Vector,
    transfer_velocity_: Vector,
) -> Vector:
    r"""
    Absolute velocity via relative and transfer velocities.

    Law:
        :code:`v_abs = v_rel + v_tr`

    Latex:
        .. math::
            {\vec v}_\text{abs} = {\vec v}_\text{rel} + {\vec v}_\text{tr}
    
    :param relative_velocity\_: velocity relative to moving frame :math:`S'`

        Symbol: :code:`v_rel`

        Latex: :math:`{\vec v}_\text{rel}`

        Dimension: *velocity*

    :param transfer_velocity\_: velocity due to movement of frame :math:`S'` relative to frame :math:`S`

        Symbol: :code:`v_tr`

        Latex: :math:`{\vec v}_\text{tr}`

        Dimension: *velocity*

    :return: velocity relative to fixed frame :math:`S`

        Symbol: :code:`v_abs`

        Latex: :math:`{\vec v}_\text{abs}`

        Dimension: *velocity*
    """

    return add_cartesian_vectors(
        relative_velocity_,
        transfer_velocity_,
    )


def relative_velocity_law(
    absolute_velocity_: Vector,
    transfer_velocity_: Vector,
) -> Vector:
    r"""
    Relative velocity via absolute and transfer velocities.

    Law:
        :code:`v_rel = v_abs - v_tr`

    Latex:
        .. math::
            {\vec v}_\text{rel} = {\vec v}_\text{abs} - {\vec v}_\text{tr}

    :param absolute_velocity\_: velocity relative to fixed frame :math:`S`

        Symbol: :code:`v_abs`

        Latex: :math:`{\vec v}_\text{abs}`

        Dimension: *velocity*
                
    :param transfer_velocity\_: velocity due to movement of frame :math:`S'` relative to frame :math:`S`

        Symbol: :code:`v_tr`

        Latex: :math:`{\vec v}_\text{tr}`

        Dimension: *velocity*
    
    :return: velocity relative to moving frame :math:`S'`

        Symbol: :code:`v_rel`

        Latex: :math:`{\vec v}_\text{rel}`

        Dimension: *velocity*
    """

    return add_cartesian_vectors(absolute_velocity_, scale_vector(-1, transfer_velocity_))


def transfer_velocity_law(
    absolute_velocity_: Vector,
    relative_velocity_: Vector,
) -> Vector:
    r"""
    Transfer velocity via absolute and relative velocities.

    Law:
        :code:`v_tr = v_abs - v_rel`

    Latex:
        .. math::
            {\vec v}_\text{tr} = {\vec v}_\text{abs} - {\vec v}_\text{rel}

    :param absolute_velocity\_: velocity relative to fixed frame :math:`S`

        Symbol: :code:`v_abs`

        Latex: :math:`{\vec v}_\text{abs}`

        Dimension: *velocity*

    :param relative_velocity\_: velocity relative to moving frame :math:`S'`

        Symbol: :code:`v_rel`

        Latex: :math:`{\vec v}_\text{rel}`

        Dimension: *velocity*
    
    :return_: velocity due to movement of frame :math:`S'` relative to frame :math:`S`

        Symbol: :code:`v_tr`

        Latex: :math:`{\vec v}_\text{tr}`

        Dimension: *velocity*
"""

    return add_cartesian_vectors(absolute_velocity_, scale_vector(-1, relative_velocity_))


@validate_input(
    relative_velocity_=units.velocity,
    transfer_velocity_=units.velocity,
)
@validate_output(units.velocity)
def calculate_absolute_velocity(
    relative_velocity_: QuantityVector,
    transfer_velocity_: QuantityVector,
) -> QuantityVector:
    result = absolute_velocity_law(
        relative_velocity_.to_base_vector(),
        transfer_velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(result)
