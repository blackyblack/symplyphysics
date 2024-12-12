"""
Photon energy is proportional to linear frequency
=================================================

Photon is the elementary part of any electromagnetical radiation, having no mass and moving with speed of light
in any reference frame. The energy of a photon depends only on its frequency.

**Notation:**

#. :quantity_notation:`planck`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Photon#Relativistic_energy_and_momentum>`__.
"""

from sympy import Eq, solve
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

frequency = symbols.temporal_frequency
"""
:symbols:`temporal_frequency` of a photon.
"""

law = Eq(energy, quantities.planck * frequency)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(photon_frequency_=frequency)
@validate_output(energy)
def calculate_energy(photon_frequency_: Quantity) -> Quantity:
    result_energy_expr = solve(law, energy, dict=True)[0][energy]
    result_expr = result_energy_expr.subs({frequency: photon_frequency_})
    return Quantity(result_expr)
