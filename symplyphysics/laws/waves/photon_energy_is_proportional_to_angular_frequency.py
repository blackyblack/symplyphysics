r"""
Photon energy is proportional to angular frequency
==================================================

Photon is the elementary part of any electromagnetical radiation, having no mass and moving with speed of light
in any reference frame. The energy of a photon depends only on its frequency.

**Notation:**

#. :math:`\hbar` is the reduced Planck constant.
"""

from sympy import Eq
from symplyphysics import (
    units,
    angle_type,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

energy = Symbol("energy", units.energy)
"""
Energy of a photon.

Symbol:
    :code:`E`
"""

angular_frequency = Symbol("angular_frequency", angle_type / units.time)
r"""
Angular angular_frequency of a photon.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

law = Eq(energy, units.hbar * angular_frequency)
r"""
:code:`E = hbar * w`

Latex:
    .. math::
        E = \hbar \omega
"""


@validate_input(angular_frequency_=angular_frequency)
@validate_output(energy)
def calculate_energy(angular_frequency_: Quantity) -> Quantity:
    result_expr = law.rhs.subs({angular_frequency: angular_frequency_})
    return Quantity(result_expr)
