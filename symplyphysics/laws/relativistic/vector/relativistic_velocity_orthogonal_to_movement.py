"""
Relativistic velocity normal to movement
========================================

Consider two inertial reference frames: one fixed (lab frame) and one tied to the moving object
(proper frame). The proper frame is moving with some velocity :math:`\\vec v` relative to the lab
frame. According to the theory of special relativity, the velocity of the object relative to lab
frame is *not* equal to the sum of its velocity in the proper frame and the velocity of the proper
frame relative to the lab frame.

**Notes:**

#. One can get the same expression for :math:`{\\vec u}_\\text{n}'` in terms of :math:`\\vec u` by
   replacing :math:`\\vec v` with :math:`-{\\vec v}`. This is essentially the inverse Lorentz
   transformation from lab frame to proper frame that uses the fact that the lab frame can be
   viewed as moving with velocity vector :math:`-{\\vec v}` relative to the proper frame.

**Conditions:**

#. Works in special relativity.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Velocity-addition_formula#General_configuration>`__.
"""

from sympy import Eq, evaluate, sqrt
from symplyphysics import validate_input, validate_output, symbols, quantities

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

normal_velocity_in_lab_frame = clone_as_vector_symbol(
    symbols.speed,
    display_symbol="u_n",
    display_latex="{\\vec u}_\\text{n}",
)
"""
Component of the velocity vector relative to the lab frame normal to :math:`\\vec v`. See
:symbols:`speed`.
"""

velocity_in_proper_frame = clone_as_vector_symbol(
    symbols.speed,
    display_symbol="u'",
    display_latex="{\\vec u'}",
)
"""
Velocity vector relative to the proper frame. See :symbols:`speed`.
"""

proper_frame_velocity = clone_as_vector_symbol(symbols.speed)
"""
Velocity vector of the proper frame relative to the lab frame. See :symbols:`speed`.
"""

with evaluate(False):
    _v = proper_frame_velocity
    _up = velocity_in_proper_frame
    _c = quantities.speed_of_light

    _up_tangential = (VectorDot(_up, _v) / VectorDot(_v, _v)) * _v
    _up_normal = _up - _up_tangential

    _numerator = sqrt(1 - VectorDot(_v, _v) / _c**2)
    _denominator = 1 + VectorDot(_up, _v) / _c**2

law = Eq(normal_velocity_in_lab_frame, _up_normal * _numerator / _denominator)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    velocity_in_proper_frame_=velocity_in_proper_frame,
    proper_frame_velocity_=proper_frame_velocity,
)
@validate_output(normal_velocity_in_lab_frame)
def calculate_orthogonal_velocity_component_in_lab_frame(
    velocity_in_proper_frame_: QuantityCoordinateVector,
    proper_frame_velocity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        velocity_in_proper_frame: velocity_in_proper_frame_,
        proper_frame_velocity: proper_frame_velocity_,
    })

    return QuantityCoordinateVector.from_expr(result)
