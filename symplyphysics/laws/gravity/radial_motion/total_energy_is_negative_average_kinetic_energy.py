"""
Total energy is negative average kinetic energy
===============================================

The total energy of an orbiting planet is equal to its negative kinetic energy averaged over time.

**Conditions:**

#. Works for elliptical (:math:`E < 0`) orbits.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

total_mechanical_energy = symbols.mechanical_energy
"""
The total :symbols:`mechanical_energy` of the planet.
"""

average_kinetic_energy = clone_as_symbol(
    symbols.kinetic_energy,
    display_symbol="avg(K)",
    display_latex="\\langle K \\rangle",
)
"""
The :symbols:`kinetic_energy` of the planet averaged over :symbols:`time`.
"""

law = Eq(total_mechanical_energy, -1 * average_kinetic_energy)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(average_kinetic_energy_=average_kinetic_energy)
@validate_output(total_mechanical_energy)
def calculate_total_energy(average_kinetic_energy_: Quantity) -> Quantity:
    result = law.rhs.subs(average_kinetic_energy, average_kinetic_energy_)
    return Quantity(result)
