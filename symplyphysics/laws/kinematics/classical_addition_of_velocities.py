r"""
Classical addition of velocities
================================

The *law of classical addition of velocities*, usually attributed to Galileo and called the
*Galilean law of velocity addition*, states that the velocity of a body in an inertial reference
frame :math:`A` can be found as a sum of its velocity in another inertial reference frame :math:`B` and the
velocity of frame :math:`B` relative to frame :math:`A`.

**Conditions:**

#. Velocity vectors must be collinear.
#. Space and time are absolute.
#. Applicable to inertial reference frames.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

body_speed_in_first_frame = clone_as_symbol(symbols.speed, subscript="OA")
"""
:symbols:`speed` of the body in frame :math:`A`.
"""

body_speed_in_second_frame = clone_as_symbol(symbols.speed, subscript="OB")
"""
:symbols:`speed` of the body in frame :math:`B`.
"""

second_frame_speed_in_first_frame = clone_as_symbol(symbols.speed, subscript="BA")
"""
:symbols:`speed` of frame :math:`B` relative to frame :math:`A`.
"""

law = Eq(
    body_speed_in_first_frame,
    body_speed_in_second_frame + second_frame_speed_in_first_frame,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    velocity_in_second_frame_=body_speed_in_second_frame,
    second_frame_velocity_in_first_frame_=second_frame_speed_in_first_frame,
)
@validate_output(body_speed_in_first_frame)
def calculate_body_velocity_in_first_frame(
    velocity_in_second_frame_: Quantity,
    second_frame_velocity_in_first_frame_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        body_speed_in_second_frame: velocity_in_second_frame_,
        second_frame_speed_in_first_frame: second_frame_velocity_in_first_frame_,
    })
    return Quantity(result)
