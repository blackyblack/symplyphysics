r"""
Photon momentum is proportional to energy
=========================================

The momentum of a photon is its energy divided by the speed of light.

See :doc:`laws.waves.photon_energy_is_proportional_to_angular_frequency` or
:doc:`laws.waves.photon_energy_is_proportional_to_linear_frequency` for energy
expressions.

**Notation:**

#. :math:`c` is the speed of light.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

momentum = Symbol("momentum", units.momentum)
"""
Momentum of a photon.

Symbol:
    :code:`p`
"""

energy = Symbol("frequency", units.energy)
"""
Energy of a photon.

Symbol:
    :code:`E`
"""

law = Eq(momentum, energy / units.speed_of_light)
r"""
:code:`p = E / c`

Latex:
    .. math::
        p = \frac{E}{c}
"""


@validate_input(photon_energy_=energy)
@validate_output(momentum)
def calculate_momentum(photon_energy_: Quantity) -> Quantity:
    result_momentum_expr = solve(law, momentum, dict=True)[0][momentum]
    result_expr = result_momentum_expr.subs({energy: photon_energy_})
    return Quantity(result_expr)
