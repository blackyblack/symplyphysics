r"""
Photon energy is proportional to linear frequency
=================================================

Photon is the elementary part of any electromagnetical radiation, having no mass and moving with speed of light
in any reference frame. The energy of a photon depends only on its frequency.

**Notation:**

#. :quantity_notation:`planck`.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, quantities)

energy = Symbol("energy", units.energy)
"""
Energy of a photon.

Symbol:
    :code:`E`
"""

frequency = Symbol("frequency", units.frequency)
r"""
Frequency of a photon.

Symbol:
    :code:`nu`

Latex:
    :math:`\nu`
"""

law = Eq(energy, quantities.planck * frequency)
r"""
:code:`E = h * nu`

Latex:
    .. math::
        E = h \nu
"""


@validate_input(photon_frequency_=frequency)
@validate_output(energy)
def calculate_energy(photon_frequency_: Quantity) -> Quantity:
    result_energy_expr = solve(law, energy, dict=True)[0][energy]
    result_expr = result_energy_expr.subs({frequency: photon_frequency_})
    return Quantity(result_expr)
