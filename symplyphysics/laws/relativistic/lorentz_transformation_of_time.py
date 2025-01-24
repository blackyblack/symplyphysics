"""
Lorentz transformation of time
==============================

Let us consider two inertial reference frames :math:`S` (lab frame) and :math:`S'` (proper frame). The proper frame moves
with speed :math:`v` relative to the lab frame. Both frames contain identical fixed (relative to their own frame)
clocks that are synchronized according to the `Einstein rule <https://en.wikipedia.org/wiki/Einstein_synchronisation>`__.
Let :math:`x, y, z, t` be the coordinates and time of some event in frame :math:`S`, and :math:`x', y', z', t'` be the coordinates
and time of the same event in frame :math:`S'`. Assuming that the space is uniform and isotropic and that the time is
uniform, there exists a linear dependence between :math:`x, y, z, t` and :math:`x', y', z', t'`, which is called the **Lorentz
transformation** of space and time.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Notes**

#. Lab frame :math:`S` is usually thought as stationary, and proper frame :math:`S'` is the one that is considered to be moving
   relative to lab frame and the moving object in question is at rest in the proper frame.
#. In this law, the Lorentz transformation from the lab frame :math:`S` into the proper frame :math:`S'` is described. In order to get
   an opposite transformation (from the proper frame :math:`S'` into the lab frame :math:`S`), replace all primed variables with unprimed
   ones and vice verce, and replace :math:`v` with :math:`-v`. This is consistent with the fact that frame :math:`S` can be viewed as moving
   with speed :math:`-v` relative to frame :math:`S'`, and hence the same Lorentz transformation can be applied.
#. In the limit :math:`v/c \\ll 1` the formula reduces to the classical Galilean transformation :math:`t' = t`.

**Conditions:**

#. Space is uniform and isotropic.
#. Time is uniform.
#. The relative frame velocity is parallel to the :math:`x`-axis.

**Links:**

#. `Wikipedia, first formula in box <https://en.wikipedia.org/wiki/Lorentz_transformation#Coordinate_transformation>`__.
"""

from sympy import Eq, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals

time_in_proper_frame = clone_as_symbol(symbols.time, display_symbol="t'", display_latex="t'")
"""
:symbols:`time` in proper frame :math:`S'`.
"""

time_in_lab_frame = symbols.time
"""
:symbols:`time` in lab frame :math:`S`.
"""

position_in_lab_frame = symbols.position
"""
:symbols:`position` in lab frame :math:`S`.
"""

proper_frame_speed_in_lab_frame = symbols.speed
"""
:symbols:`speed` of proper frame :math:`S'` relative to lab frame :math:`S`.
"""

law = Eq(time_in_proper_frame, (time_in_lab_frame -
    proper_frame_speed_in_lab_frame * position_in_lab_frame / quantities.speed_of_light**2) / sqrt(1 -
    (proper_frame_speed_in_lab_frame / quantities.speed_of_light)**2))
"""
:laws:symbol::

:laws:latex::
"""

# In the limit `V/c << 1`, the formula reduces to the classical relation `t' = t`.
_classical_time_expr = law.rhs.series(proper_frame_speed_in_lab_frame, 0, 1).removeO()
assert expr_equals(_classical_time_expr, time_in_lab_frame)


@validate_input(
    time_in_lab_frame_=time_in_lab_frame,
    position_in_lab_frame_=position_in_lab_frame,
    proper_frame_speed_in_lab_frame_=proper_frame_speed_in_lab_frame,
)
@validate_output(time_in_proper_frame)
def calculate_time_in_proper_frame(
    time_in_lab_frame_: Quantity,
    position_in_lab_frame_: Quantity,
    proper_frame_speed_in_lab_frame_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        time_in_lab_frame: time_in_lab_frame_,
        position_in_lab_frame: position_in_lab_frame_,
        proper_frame_speed_in_lab_frame: proper_frame_speed_in_lab_frame_,
    })
    return Quantity(result)
