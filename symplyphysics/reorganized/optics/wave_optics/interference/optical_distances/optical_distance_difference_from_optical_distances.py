"""
Optical distance difference from optical distances
==================================================

The optical difference in the course of two rays is the difference
in the optical distances traversed by each of the rays

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Optical_path_length#Optical_path_difference>`__.

..
    TODO elaborate on optical distance
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

first_optical_distance = clone_as_symbol(symbols.optical_distance, subscript="1")
"""
:symbols:`optical_distance` of first wave.
"""

second_optical_distance = clone_as_symbol(symbols.optical_distance, subscript="2")
"""
:symbols:`optical_distance` of second wave.
"""

optical_difference_distance = clone_as_symbol(
    symbols.optical_distance,
    display_symbol="Delta(Lambda)",
    display_latex="\\Delta \\Lambda",
)
"""
:symbols:`optical_distance` difference between waves.
"""

law = Eq(optical_difference_distance, second_optical_distance - first_optical_distance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(optical_distance1_=first_optical_distance,
    optical_distance2_=second_optical_distance)
@validate_output(optical_difference_distance)
def calculate_optical_difference_distance(optical_distance1_: Quantity,
    optical_distance2_: Quantity) -> Quantity:
    solved = solve(law, optical_difference_distance, dict=True)[0][optical_difference_distance]
    result_expr = solved.subs({
        first_optical_distance: optical_distance1_,
        second_optical_distance: optical_distance2_
    })
    return Quantity(result_expr)
