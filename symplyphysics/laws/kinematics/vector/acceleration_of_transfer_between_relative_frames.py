"""
Acceleration of transfer between relative frames
================================================

Imagine two reference frames, one of which is fixed (:math:`S`) and the other one is moving
(:math:`S'`). The motion of a body stationary in moving frame :math:`S'` due to the motion of the
frame itself is called transfer motion. The acceleration related to such motion is called transfer
acceleration. It is composed of the acceleration of the moving frame relative to the fixed frame,
centripetal acceleration and the acceleration due to uneven rotation of the moving frame. The
transfer acceleration only depends on the motion of frame :math:`S'` relative to stationary frame
:math:`S`, so its physical meaning would be that it is the acceleration in :math:`S` of a point
stationary in :math:`S'`.

**Links:**

#. `Wikipedia <https://ru.wikipedia.org/wiki/%D0%A1%D0%BB%D0%BE%D0%B6%D0%BD%D0%BE%D0%B5_%D0%B4%D0%B2%D0%B8%D0%B6%D0%B5%D0%BD%D0%B8%D0%B5#%D0%A3%D1%81%D0%BA%D0%BE%D1%80%D0%B5%D0%BD%D0%B8%D0%B5>`__.

..
    TODO find English link
"""

from sympy import Eq
from symplyphysics import validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

transfer_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_tr",
    display_latex="{\\vec a}_\\text{tr}",
)
"""
Vector of the transfer :symbols:`acceleration`.
"""

moving_frame_acceleration = clone_as_vector_symbol(symbols.acceleration, subscript="0")
"""
Vector of the :symbols:`acceleration` of the :math:`S'` relative to :math:`S`.
"""

centripetal_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_cp",
    display_latex="{\\vec a}_\\text{cp}",
)
"""
Vector of the body's centripetal :symbols:`acceleration` relative to :math:`S'`.
"""

rotation_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_rot",
    display_latex="{\\vec a}_\\text{rot}",
)
"""
Vector of acceleration due to non-uniform rotation of :math:`S'`.
"""

law = Eq(
    transfer_acceleration,
    moving_frame_acceleration + centripetal_acceleration + rotation_acceleration,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    moving_frame_acceleration_=moving_frame_acceleration,
    centripetal_acceleration_=centripetal_acceleration,
    rotation_acceleration_=rotation_acceleration,
)
@validate_output(transfer_acceleration)
def calculate_transfer_acceleration(
    moving_frame_acceleration_: QuantityCoordinateVector,
    centripetal_acceleration_: QuantityCoordinateVector,
    rotation_acceleration_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        moving_frame_acceleration: moving_frame_acceleration_,
        centripetal_acceleration: centripetal_acceleration_,
        rotation_acceleration: rotation_acceleration_,
    })

    return QuantityCoordinateVector.from_expr(result)
