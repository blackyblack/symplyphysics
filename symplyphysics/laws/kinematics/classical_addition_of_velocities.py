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
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

body_speed_in_first_frame = Symbol("body_speed_in_first_frame", units.velocity)
"""
Speed of the body in frame :math:`A`.

Symbol:
    :code:`u`
"""

body_speed_in_second_frame = Symbol("body_speed_in_second_frame", units.velocity)
r"""
Speed of the body in frame :math:`B`.

Symbol:
    :code:`u'`

Latex:
    :math:`u'`
"""

second_frame_speed_in_first_frame = Symbol("second_frame_speed_in_first_frame", units.velocity)
"""
Speed of frame :math:`B` relative to frame :math:`A`.

Symbol:
    :code:`v`
"""

law = Eq(
    body_speed_in_first_frame,
    body_speed_in_second_frame + second_frame_speed_in_first_frame,
)
r"""
:code:`u = u' + v`

Latex:
    .. math::
        u = u' + v
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
