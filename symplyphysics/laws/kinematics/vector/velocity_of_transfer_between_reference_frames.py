"""
Velocity of transfer between reference frames
=============================================

Suppose two reference frames, one of which is fixed (:math:`S`) and the other one is moving
(:math:`S'`). The movement of a body stationary in moving frame :math:`S'` due to the movement of
the frame itself is called transfer movement. The velocity related to such movement is called
transfer velocity. For any material point :math:`X`, its transfer velocity relative to fixed frame
:math:`S` is the sum of the velocity of frame :math:`S'` relative to frame :math:`S` and the cross
product of the angular velocity of moving frame's rotation and the position vector of :math:`X` in
moving frame :math:`S'`.

**Links:**

#. `Wikipedia, first formula <https://ru.wikipedia.org/wiki/%D0%A2%D0%B5%D0%BE%D1%80%D0%B5%D0%BC%D0%B0_%D0%BE_%D1%81%D0%BB%D0%BE%D0%B6%D0%B5%D0%BD%D0%B8%D0%B8_%D1%81%D0%BA%D0%BE%D1%80%D0%BE%D1%81%D1%82%D0%B5%D0%B9#%D0%9E%D0%B1%D1%81%D1%83%D0%B6%D0%B4%D0%B5%D0%BD%D0%B8%D0%B5>`__.

..
    TODO: find English link
"""

from sympy import Eq
from symplyphysics import validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorCross
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.solvers import solve_for_vector

transfer_velocity = clone_as_vector_symbol(
    symbols.speed,
    display_symbol="v_tr",
    display_latex="{\\vec v}_\\text{tr}",
)
"""
Vector of transfer velocity of point :math:`X` relative to fixed frame :math:`S`. See :symbols:`speed`.
"""

moving_frame_velocity = clone_as_vector_symbol(symbols.speed, subscript="0")
"""
Vector of moving frame :math:`S'` relative to fixed frame :math:`S`.
"""

angular_velocity = clone_as_vector_symbol(symbols.angular_speed)
"""
Pseudovector of the angular velocity related to the rotation of moving frame :math:`S'` about the
instantaneous axis. See :symbols:`angular_speed`.
"""

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Position vector of point :math:`X` relative to moving frame :math:`S'`. See
:symbols:`distance_to_origin`.
"""

law = Eq(
    transfer_velocity,
    moving_frame_velocity + VectorCross(angular_velocity, position_vector),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    moving_frame_velocity_=moving_frame_velocity,
    angular_velocity_=angular_velocity,
    position_vector_=position_vector,
)
@validate_output(transfer_velocity)
def calculate_transfer_velocity(
    moving_frame_velocity_: QuantityCoordinateVector,
    angular_velocity_: QuantityCoordinateVector,
    position_vector_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = solve_for_vector(law, transfer_velocity).subs({
        moving_frame_velocity: moving_frame_velocity_,
        angular_velocity: angular_velocity_,
        position_vector: position_vector_,
    })

    return QuantityCoordinateVector.from_expr(result)
