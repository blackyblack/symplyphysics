"""
Relativistic length via proper length and speed
===============================================

Length contraction is the phenomenon that a moving object's length is measured to be shorter
than its proper length, which is the length as measured in the object's own rest frame. Note
that this phenomenon is only observed in the direction parallel to the velocity of the object.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Length_contraction#Basis_in_relativity>`__.

..
    TODO rename file
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities,
    clone_as_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.special_relativity.relativistic_kinematics.lorentz_transformation import lorentz_transformation_of_coordinate as coordinate_transformation

proper_length = symbols.proper_length
"""
:symbols:`proper_length` of the object.
"""

speed = symbols.speed
"""
:symbols:`speed` of the object.
"""

relativistic_length = symbols.length
"""
Relativistic :symbols:`length` of the object, i.e. its length measured in the external reference
frame.
"""

law = Eq(relativistic_length, proper_length * sqrt(1 - speed**2 / quantities.speed_of_light**2))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law from the Lorentz transformation of coordinate.

## Consider a rod at rest in the proper frame S', which moves with the given speed relative
## to the lab frame. The length of the rod in the lab frame is measured by locating both of
## its ends at the same lab time, so the positions of the ends in the lab frame differ by
## the relativistic length.
_first_end_position = clone_as_symbol(symbols.position, subscript="1")
_second_end_position = _first_end_position + relativistic_length

## Transform the positions of both ends into the proper frame S'. Since both measurement
## events happen at the same lab time, the time-dependent terms cancel out in the difference.
_first_end_position_transformed = coordinate_transformation.law.rhs.subs(
    coordinate_transformation.position_in_lab_frame, _first_end_position)
_second_end_position_transformed = coordinate_transformation.law.rhs.subs(
    coordinate_transformation.position_in_lab_frame, _second_end_position)

## The distance between the ends of the rod in its rest frame S' is the proper length.
_proper_length_expr = (_second_end_position_transformed - _first_end_position_transformed).subs(
    coordinate_transformation.proper_frame_speed_in_lab_frame, speed)

_relativistic_length_derived = solve(Eq(proper_length, _proper_length_expr),
    relativistic_length)[0]

assert expr_equals(_relativistic_length_derived, law.rhs)


@validate_input(rest_length_=proper_length, velocity_=speed)
@validate_output(relativistic_length)
def calculate_relativistic_length(rest_length_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, relativistic_length, dict=True)[0][relativistic_length]
    length_applied = result_expr.subs({proper_length: rest_length_, speed: velocity_})
    return Quantity(length_applied)
