"""
Photon energy is proportional to angular frequency
==================================================

Photon is the elementary part of any electromagnetical radiation, having no mass and moving with speed of light
in any reference frame. The energy of a photon depends only on its frequency.

**Notation:**

#. :quantity_notation:`hbar`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Photon#Relativistic_energy_and_momentum>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    quantities,
    symbols,
)

energy = symbols.energy
"""
:symbols:`energy` of a photon.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of a photon.
"""

law = Eq(energy, quantities.hbar * angular_frequency)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(angular_frequency_=angular_frequency)
@validate_output(energy)
def calculate_energy(angular_frequency_: Quantity) -> Quantity:
    result_expr = law.rhs.subs({angular_frequency: angular_frequency_})
    return Quantity(result_expr)
