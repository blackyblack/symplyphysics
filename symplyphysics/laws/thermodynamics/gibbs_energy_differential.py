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
"""

from sympy import Eq
from symplyphysics import (
    clone_as_symbol,
    symbols,
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

gibbs_energy_change = Symbol("gibbs_energy_change", units.energy)
"""
Infinitesimal change in Gibbs energy of the system.

Symbol:
    :code:`dG`
"""

entropy = Symbol("entropy", units.energy / units.temperature)
"""
Entropy of the system.

Symbol:
    :code:`S`
"""

temperature_change = clone_as_symbol(symbols.temperature, display_symbol="dT")
"""
Infinitesimal change in :symbols:`temperature` of the system.
"""

volume = Symbol("volume", units.volume)
"""
Volume of the system.

Symbol:
    :code:`V`
"""

pressure_change = Symbol("pressure_change", units.pressure)
"""
Infinitesimal change in pressure inside the system.

Symbol:
    :code:`dp`
"""

chemical_potential = Symbol("chemical_potential", units.energy)
r"""
Chemical potential of the system.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

particle_count_change = Symbol("particle_count_change", dimensionless)
"""
Infinitesimal change in the number of particles in the system.

Symbol:
    :code:`dN`
"""

law = Eq(
    gibbs_energy_change,
    -1 * entropy * temperature_change + volume * pressure_change +
    chemical_potential * particle_count_change,
)
r"""
:code:`dG = -1 * S * dT + V * dp + mu * dN`

Latex:
    .. math::
        dG = - S \, dT + V \, dp + \mu \, dN
"""


#pylint: disable=too-many-arguments
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
    result = law.rhs.subs({
        entropy: entropy_,
        temperature_change: temperature_change_,
        volume: volume_,
        pressure_change: pressure_change_,
        chemical_potential: chemical_potential_,
        particle_count_change: particle_count_change_,
    })
    return Quantity(result)
