"""
Velocity relative to reference frame
====================================

For any reference frame, whether it is inertial or not, the motion relative to it can be described using
the position vector relative to that frame's origin.
"""

from sympy import Expr, symbols
from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    Quantity,
    QuantityVector,
    add_cartesian_vectors,
    subtract_cartesian_vectors,
    scale_vector,
)
from symplyphysics.core.vectors.arithmetics import (
    diff_cartesian_vector,
    integrate_cartesian_vector,
)

def relative_velocity_law(
    position_: Vector,
    time_: Expr,
) -> Vector:
    r"""
    Velocity relative to :math:`S`.

    Law:
        :code:`v_rel = Derivative(r(t), t)`

    Latex:
        .. math::
            {\vec v}_\text{rel} = \frac{d \vec r}{d t}

    :param position\_: radius vector, or position vector, of body in :math:`S` as a function of time

        Symbol: :code:`r(t)`

        Latex: :math:`\vec r(t)`

        Dimension: *length*

    :param time\_: time

        Symbol: :code:`t`

        Dimension: *time*

    :return: velocity relative to :math:`S`

        Symbol: :code:`v_rel`

        Latex: :math:`{\vec v}_\text{rel}`

        Dimension: *velocity*
    """

    return diff_cartesian_vector(position_, time_)


def relative_position_law(
    initial_position_: Vector,
    velocity_: Vector,
    time_: Expr,
) -> Vector:
    r"""
    Final position via initial position and velocity as a function of time.

    Law:
        :code:`r = r_0 + Integral(v_rel(t), t)`

    Latex:
        .. math::
            \vec r = {\vec r}_0 + \int {\vec v}_\text{rel}(t) dt

    :param initial\_position\_: position vector in :math:`S` at :math:`t = 0`

        Symbol: :code:`r_0`

        Latex: :math:`{\vec r}_0`

        Dimension: *length*

    :param velocity\_: velocity relative to :math:`S` as a function of time

        Symbol: :code:`v_rel(t)`

        Latex: :math:`{\vec v}_\text{rel}(t)`

        Dimension: *velocity*

    :param time\_: time

        Symbol: :code:`t`

        Dimension: *time*

    :return: position vector in :math:`S` at time :math:`t`

        Symbol: :code:`r`

        Latex: :math:`\vec r`

        Dimension: *length*
    """

    return add_cartesian_vectors(
        initial_position_,
        integrate_cartesian_vector(velocity_, time_),
    )


@validate_input(
    position_before_=units.length,
    position_after_=units.length,
    time_change_=units.time,
)
@validate_output(units.velocity)
def calculate_relative_velocity(
    position_before_: QuantityVector,
    position_after_: QuantityVector,
    time_change_: Quantity,
) -> QuantityVector:
    time_ = symbols("time")
    position_ = scale_vector(
        time_ / time_change_,
        subtract_cartesian_vectors(
        position_after_.to_base_vector(),
        position_before_.to_base_vector(),
        ),
    )
    velocity_ = relative_velocity_law(position_, time_)
    return QuantityVector.from_base_vector(velocity_)
