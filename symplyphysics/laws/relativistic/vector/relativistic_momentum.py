"""
Relativistic momentum via rest mass and velocity
================================================

Momentum (amount of motion) is a vector physical quantity that is a measure of the mechanical
movement of a body. The relativistic momentum also takes into account speed limits equal to the
speed of light.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Notes:**

#. To find rest mass refer to the :ref:`scalar law <Relativistic momentum via rest mass and
   speed>`.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Mass_in_special_relativity#Relativistic_mass>`__.
"""

from sympy import Eq, sqrt
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.quantities import speed_of_light

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

momentum = clone_as_vector_symbol(symbols.momentum)
"""
Vector of the body's relativistic :symbols:`momentum`.
"""

rest_mass = symbols.rest_mass
"""
:symbols:`rest_mass` of the body.
"""

velocity = clone_as_vector_symbol(symbols.speed)
"""
Vector of the body's velocity. See :symbols:`speed`.
"""

momentum_law_ = Eq(
    momentum,
    (rest_mass * velocity) / sqrt(1 - VectorDot(velocity, velocity) / speed_of_light**2),
)
"""
:laws:symbol::

:laws:latex::
"""

velocity_law_ = Eq(
    velocity,
    (momentum * speed_of_light) /
    sqrt((rest_mass * speed_of_light)**2 + VectorDot(momentum, momentum)),
)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: prove that the expressions for momentum and velocity are equivalent.


@validate_input(
    rest_mass_=rest_mass,
    velocity_=velocity,
)
@validate_output(momentum)
def calculate_momentum(
    rest_mass_: Quantity,
    velocity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = momentum_law_.rhs.subs({
        rest_mass: rest_mass_,
        velocity: velocity_,
    })

    return QuantityCoordinateVector.from_expr(result)
