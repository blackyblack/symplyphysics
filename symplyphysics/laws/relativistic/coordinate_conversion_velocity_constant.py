"""
Coordinate conversion at constant velocity
==========================================

Let the frame of reference move relative to the other frame of reference at a constant speed along the X
axis, and the origin of the spatial coordinates coincide at the initial moment of time in both systems.
Then there are simple transformations that can be used to get the x coordinate in one frame of reference,
knowing the x coordinate in another frame of reference.

**Links:**

#. `Wikipedia, second formula in box <https://en.wikipedia.org/wiki/Lorentz_transformation#Coordinate_transformation>`__.

..
    TODO rename file
"""

from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    quantities,
    symbols,
    clone_as_symbol,
)

position_in_second_frame = clone_as_symbol(symbols.position, subscript="2")
"""
:symbols:`position` in the second frame of reference.
"""

position_in_first_frame = clone_as_symbol(symbols.position, subscript="1")
"""
:symbols:`position` in the first frame of reference.
"""

second_frame_speed_in_first_frame = symbols.speed
"""
:symbols:`speed` of the second reference frame relative to the first one.
"""

time_in_first_frame = clone_as_symbol(symbols.time, subscript="1")
"""
:symbols:`time` in the first frame of reference.
"""

law = Eq(position_in_second_frame, (position_in_first_frame - second_frame_speed_in_first_frame * time_in_first_frame) / sqrt(
    (1 - (second_frame_speed_in_first_frame / quantities.speed_of_light)**2)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(coordinate_first_frame_=position_in_first_frame,
    velocity_=second_frame_speed_in_first_frame,
    time_first_frame_=time_in_first_frame)
@validate_output(position_in_second_frame)
def calculate_coordinate_second_frame(coordinate_first_frame_: Quantity, velocity_: Quantity,
    time_first_frame_: Quantity) -> Quantity:
    result_coordinate_first_frame_second_frame_expr = solve(law, position_in_second_frame,
        dict=True)[0][position_in_second_frame]
    result_expr = result_coordinate_first_frame_second_frame_expr.subs({
        position_in_first_frame: coordinate_first_frame_,
        second_frame_speed_in_first_frame: velocity_,
        time_in_first_frame: time_first_frame_
    })
    return Quantity(result_expr)
