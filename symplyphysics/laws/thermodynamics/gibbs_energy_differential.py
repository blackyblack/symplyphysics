r"""
Gibbs energy differential
=========================

The fundamental thermodynamic relations are fundamental equations which demonstate how important
thermodynamic quantities depend on variables that are measurable experimentally.

**Notation:**

#. :math:`d` denotes an exact, path-independent differential.

**Notes:**

#. Temperature, pressure, and particle count are so called natural variables of Gibbs energy as a
   thermodynamic potential.
#. For a system with more than one type of particles, the last term can be represented as a sum over all
   types of particles, i.e. :math:`\sum_i \mu_i \, d N_i`.

**Conditions:**

#. The system is in thermal equilibrium with its surroundings.
#. The system is composed of only one type of particles, i.e. the system is a pure substance.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Fundamental_thermodynamic_relation>`__.
"""

from sympy import Eq
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

gibbs_energy_change = clone_as_symbol(symbols.gibbs_energy, display_symbol="dG", display_latex="dG")
"""
Infinitesimal change in :symbols:`gibbs_energy` of the system.
"""

entropy = symbols.entropy
"""
:symbols:`entropy` of the system.
"""

temperature_change = clone_as_symbol(symbols.temperature, display_symbol="dT", display_latex="dT")
"""
Infinitesimal change in :symbols:`temperature` of the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

pressure_change = clone_as_symbol(symbols.pressure, display_symbol="dp", display_latex="dp")
"""
Infinitesimal change in :symbols:`pressure` inside the system.
"""

chemical_potential = symbols.chemical_potential
r"""
:symbols:`chemical_potential` of the system.
"""

particle_count_change = clone_as_symbol(symbols.particle_count, display_symbol="dN", display_latex="dN")
"""
Infinitesimal change in the :symbols:`particle_count` of the system.
"""

law = Eq(
    gibbs_energy_change,
    -1 * entropy * temperature_change + volume * pressure_change +
    chemical_potential * particle_count_change,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    entropy_=entropy,
    temperature_change_=temperature_change,
    volume_=volume,
    pressure_change_=pressure_change,
    chemical_potential_=chemical_potential,
    particle_count_change_=particle_count_change,
)
@validate_output(gibbs_energy_change)
def calculate_gibbs_energy_change(
    entropy_: Quantity,
    temperature_change_: Quantity,
    volume_: Quantity,
    pressure_change_: Quantity,
    chemical_potential_: Quantity,
    particle_count_change_: int,
) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    result = law.rhs.subs({
        entropy: entropy_,
        temperature_change: temperature_change_,
        volume: volume_,
        pressure_change: pressure_change_,
        chemical_potential: chemical_potential_,
        particle_count_change: particle_count_change_,
    })
    return Quantity(result)
