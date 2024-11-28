"""
Average potential energy via average kinetic energy
===================================================

The average potential energy of an orbiting planet is directly proportional to its average
kinetic energy. The average of both energies is taken w.r.t. time.

**Conditions:**

#. Works for elliptical orbits, in which case the total energy :math:`E` of the planet is negative.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

average_potential_energy = clone_as_symbol(
    symbols.potential_energy,
    display_symbol="avg(U)",
    display_latex="\\langle U \\rangle",
)
"""
The :symbols:`potential_energy` of the planet averaged over :symbols:`time`.
"""

average_kinetic_energy = clone_as_symbol(
    symbols.kinetic_energy,
    display_symbol="avg(K)",
    display_latex="\\langle K \\rangle",
)
"""
The :symbols:`kinetic_energy` of the planet averaged over :symbols:`time`.
"""

law = Eq(average_potential_energy, -2 * average_kinetic_energy)
"""
:laws:symbol::

:laws:latex::
"""

# TODO Prove using the virial's theorem.


@validate_input(average_kinetic_energy_=average_kinetic_energy)
@validate_output(average_potential_energy)
def calculate_average_potential_energy(average_kinetic_energy_: Quantity) -> Quantity:
    result = law.rhs.subs(average_kinetic_energy, average_kinetic_energy_)
    return Quantity(result)
