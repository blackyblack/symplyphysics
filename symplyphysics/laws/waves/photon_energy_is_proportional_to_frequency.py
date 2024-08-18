"""
Photon energy is proportional to linear frequency
=================================================

Photon is the elementary part of any electromagnetical radiation, having no mass and moving with speed of light
in any reference frame. The energy of a photon depends only on its frequency.
"""

from sympy import (Eq, solve)
from sympy.physics.units import planck as planck_constant
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

photon_energy = Symbol("photon_energy", units.energy)
r"""
Energy of a photon.

Symbol:
    :code:`E_ph`

Latex:
    :math:`E_\text{ph}`
"""

photon_frequency = Symbol("frequency", units.frequency)
r"""
Frequency of a photon.

Symbol:
    :code:`nu`

Latex:
    :math:`\nu`
"""

law = Eq(photon_energy, planck_constant * photon_frequency)
r"""
:code:`E_ph = h * nu`

Latex:
    .. math::
        E_\text{ph} = h \nu
"""


@validate_input(photon_frequency_=photon_frequency)
@validate_output(photon_energy)
def calculate_energy(photon_frequency_: Quantity) -> Quantity:
    result_energy_expr = solve(law, photon_energy, dict=True)[0][photon_energy]
    result_expr = result_energy_expr.subs({photon_frequency: photon_frequency_})
    return Quantity(result_expr)
