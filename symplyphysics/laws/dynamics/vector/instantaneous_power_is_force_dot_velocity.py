"""
Instantaneous power is dot product of force and velocity
========================================================

**Instantaneous power** is a measure of the rate of energy transfer or conversion at a given point
in time.

**Conditions:**

#. Force :math:`\\vec F` is constant.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Power_(physics)#Definition>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

power = symbols.power
"""
Instantaneous :symbols:`power` due to :attr:`~force`.
"""

force = clone_as_vector_symbol(symbols.force)
"""
Vector of :symbols:`force` exerted on the body.
"""

velocity = clone_as_vector_symbol(symbols.speed)
"""
Vector of the body's velocity. See :symbols:`speed`.
"""

law = Eq(power, VectorDot(force, velocity))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(force_=force, velocity_=velocity)
@validate_output(power)
def calculate_power(
    force_: QuantityCoordinateVector,
    velocity_: QuantityCoordinateVector,
) -> Quantity:
    result = law.rhs.subs({
        force: force_,
        velocity: velocity_,
    })
    return Quantity(result)
