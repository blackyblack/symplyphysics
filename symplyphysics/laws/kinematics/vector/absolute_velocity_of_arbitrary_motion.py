"""
Absolute velocity of arbitrary motion
=====================================

Imagine two reference frames, one of which is fixed (:math:`S`) and the other one is moving
arbitrarily (:math:`S'`). The motion of the body relative to fixed frame :math:`S` is called
*absolute motion*. The motion of the body relative to moving frame :math:`S'` is called *relative
motion*. The motion of the body due to the motion of reference frame :math:`S'` is called *transfer
motion*. Absolute velocity is the sum of relative and transfer velocities.

**Notes:**

#. Moving frame :math:`S'` can perform both translational and rotational motion.

**Links:**

#. `Wikipedia <https://ru.wikipedia.org/wiki/%D0%A1%D0%BB%D0%BE%D0%B6%D0%BD%D0%BE%D0%B5_%D0%B4%D0%B2%D0%B8%D0%B6%D0%B5%D0%BD%D0%B8%D0%B5#%D0%A1%D0%BA%D0%BE%D1%80%D0%BE%D1%81%D1%82%D1%8C>`__.

..
    TODO: find English link
"""

from sympy import Eq
from symplyphysics import validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

absolute_velocity = clone_as_vector_symbol(
    symbols.speed,
    display_symbol="v_abs",
    display_latex="{\\vec v}_\\text{abs}",
)
"""
Vector of velocity relative to fixed frame :math:`S`. See :symbols:`speed`.
"""

relative_velocity = clone_as_vector_symbol(
    symbols.speed,
    display_symbol="v_rel",
    display_latex="{\\vec v}_\\text{rel}",
)
"""
Vector of velocity relative to moving frame :math:`S'`. See :symbols:`speed`.
"""

transfer_velocity = clone_as_vector_symbol(
    symbols.speed,
    display_symbol="v_tr",
    display_latex="{\\vec v}_\\text{tr}",
)
"""
Vector of velocity due to movement of :math:`S'` relative to :math:`S`. See :symbols:`speed`.
"""

law = Eq(absolute_velocity, relative_velocity + transfer_velocity)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    relative_velocity_=relative_velocity,
    transfer_velocity_=transfer_velocity,
)
@validate_output(absolute_velocity)
def calculate_absolute_velocity(
    relative_velocity_: QuantityCoordinateVector,
    transfer_velocity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = law.rhs.subs({
        relative_velocity: relative_velocity_,
        transfer_velocity: transfer_velocity_,
    })

    return QuantityCoordinateVector.from_expr(result)
