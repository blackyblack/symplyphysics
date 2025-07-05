"""
Coriolis acceleration
=====================

Imagine two reference frames, one of which is fixed (:math:`S`) and the other one is moving
(:math:`S'`). When the body is moving within a rotating coordinate system, its path deflects due
to the appearance of the Coriolis acceleration on it. The object does not actually deviate from
its path per se but it appears to do so because of the motion of the coordinate system.

Suppose a reference frame :math:`S'` is fixed to a rotating body :math:`A` (e.g. Earth), so that
frame :math:`S'` rotates w.r.t. another static reference frame :math:`S`. The Coriolis
acceleration is the acceleration another body :math:`B` has when moving within rotating reference
frame :math:`S'`, so it is essentially zero for objects at rest in :math:`S'`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Coriolis_force#Formula>`__.

..
    TODO: rename law to include its dependency on relative and angular velocity
"""

from sympy import Eq
from symplyphysics import validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorCross
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

coriolis_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_Cor",
    display_latex="{\\vec a}_\\text{Cor}",
)
"""
Vector of the body's Coriolis :symbols:`acceleration` in :math:`S'`.
"""

relative_velocity = clone_as_vector_symbol(
    symbols.speed,
    display_symbol="v_rel",
    display_latex="{\\vec v}_\\text{rel}",
)
"""
Vector of the body's velocity relative to :math:`S'`. See :symbols:`speed`.
"""

angular_velocity = clone_as_vector_symbol(symbols.angular_speed)
"""
Pseudovector of the angular velocity of the body's rotation. See :symbols:`angular_speed`.
"""

law = Eq(coriolis_acceleration, 2 * VectorCross(relative_velocity, angular_velocity))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    angular_velocity_=angular_velocity,
    velocity_=relative_velocity,
)
@validate_output(coriolis_acceleration)
def calculate_coriolis_acceleration(
    angular_velocity_: QuantityCoordinateVector,
    velocity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        relative_velocity: velocity_,
        angular_velocity: angular_velocity_,
    })

    return QuantityCoordinateVector.from_expr(result)
