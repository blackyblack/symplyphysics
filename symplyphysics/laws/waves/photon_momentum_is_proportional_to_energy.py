"""
Photon momentum is proportional to energy
=========================================

The momentum of a photon is its energy divided by the speed of light.

See :doc:`laws.waves.photon_energy_is_proportional_to_angular_frequency` or
:doc:`laws.waves.photon_energy_is_proportional_to_linear_frequency` for energy
expressions.

**Notation:**

#. :quantity_notation:`speed_of_light`.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, quantities, symbols

momentum = symbols.momentum
"""
:symbols:`momentum` of a photon.
"""

energy = symbols.energy
"""
:symbols:`energy` of a photon.
"""

law = Eq(momentum, energy / quantities.speed_of_light)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(photon_energy_=energy)
@validate_output(momentum)
def calculate_momentum(photon_energy_: Quantity) -> Quantity:
    result_momentum_expr = solve(law, momentum, dict=True)[0][momentum]
    result_expr = result_momentum_expr.subs({energy: photon_energy_})
    return Quantity(result_expr)
