r"""
Kinetic energy via angular momentum and angular velocity
========================================================

Kinetic energy of a body rotating around a fixed or instantaneous axis depends on its
angular momentum and angular velocity.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

kinetic_energy = symbols.kinetic_energy
"""
:symbols:`kinetic_energy` of the rotating body.
"""

angular_momentum = clone_as_vector_symbol(symbols.angular_momentum)
"""
Pseudovector of the body's :symbols:`angular_momentum`.
"""

angular_velocity = clone_as_vector_symbol(symbols.angular_speed)
"""
Pseudovector of the body's angular velocity. See :symbols:`angular_speed`.
"""

law = Eq(kinetic_energy, VectorDot(angular_momentum, angular_velocity) / 2)
"""
:laws:symbols::

:laws:latex::
"""


@validate_input(
    angular_momentum_=angular_momentum,
    angular_velocity_=angular_velocity,
)
@validate_output(kinetic_energy)
def calculate_kinetic_energy(
    angular_momentum_: QuantityCoordinateVector,
    angular_velocity_: QuantityCoordinateVector,
) -> Quantity:
    result = law.rhs.subs({
        angular_momentum: angular_momentum_,
        angular_velocity: angular_velocity_,
    })

    return Quantity(result)
