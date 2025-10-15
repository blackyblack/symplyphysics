"""
Energy—momentum relation
========================

Relativistic momentum and total relativistic energy of a relativistic particle are related by a
linear equation.

**Conditions:**

#. Velocity :math:`\\vec v` and momentum :math:`\\vec p` must be parallel to each other.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Energy%E2%80%93momentum_relation#Heuristic_approach_for_massive_particles>`__.

..
    TODO: find a more exact link
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.quantities import speed_of_light

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorNorm
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.solvers import solve_for_vector

total_energy = symbols.energy
"""
Total energy of the relativistic :symbols:`energy`.
"""

momentum = clone_as_vector_symbol(symbols.momentum)
"""
Vector of the particle's relativistic :symbols:`momentum`.
"""

velocity = clone_as_vector_symbol(symbols.speed)
"""
Vector of the particle's velocity. See :symbols:`speed`.
"""

vector_law = Eq(momentum * speed_of_light**2, total_energy * velocity)
"""
:laws:symbol::

:laws:latex::
"""

energy_law = Eq(
    total_energy,
    speed_of_light**2 * (VectorNorm(momentum) / VectorNorm(velocity)),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    total_energy_=total_energy,
    velocity_=velocity,
)
@validate_output(momentum)
def calculate_momentum(
    total_energy_: Quantity,
    velocity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = solve_for_vector(vector_law, momentum).subs({
        total_energy: total_energy_,
        velocity: velocity_,
    })

    return QuantityCoordinateVector.from_expr(result)


# UNIQUE_LAW_ID: 711
