"""
Relativistic sum of velocities
==============================

In relativistic mechanics, the speed of the body in the lab reference frame
is no longer the sum of the speed of the body's proper reference frame and
its speed in the proper frame.

**Conditions:**

#. The velocity of the body relative to the proper frame and the velocity of the
   proper frame relative to the lab frame must be parallel to each other.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Velocity-addition_formula#Special_relativity>`__.
"""

from sympy import Eq, solve
from symplyphysics import (Quantity, validate_input,
    validate_output, symbols, quantities, clone_as_symbol)

body_speed_in_lab_frame = clone_as_symbol(symbols.speed, subscript="OL")
"""
:symbols:`speed` of the body relative to the lab frame.
"""

body_speed_in_proper_frame = clone_as_symbol(symbols.speed, subscript="OP")
"""
:symbols:`speed` of the body relative to the proper reference frame.
"""

proper_frame_speed_in_lab_frame = clone_as_symbol(symbols.speed, subscript="PL")
"""
:symbols:`speed` of the proper frame relative to the lab frame.
"""

law = Eq(body_speed_in_lab_frame, (body_speed_in_proper_frame + proper_frame_speed_in_lab_frame) / (1 +
    (body_speed_in_proper_frame * proper_frame_speed_in_lab_frame) / quantities.speed_of_light**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    first_velocity_=body_speed_in_proper_frame,
    second_velocity_=proper_frame_speed_in_lab_frame,
)
@validate_output(body_speed_in_lab_frame)
def calculate_velocity(first_velocity_: Quantity, second_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, body_speed_in_lab_frame)[0]
    velocity_applied = result_expr.subs({
        body_speed_in_proper_frame: first_velocity_,
        proper_frame_speed_in_lab_frame: second_velocity_,
    })
    return Quantity(velocity_applied)
