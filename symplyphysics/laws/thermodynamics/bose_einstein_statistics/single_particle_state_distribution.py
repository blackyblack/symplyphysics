r"""
Single particle state distribution
==================================

Occupancy of a single-particle state of bosons is the probability of a single boson
to occupy a state with a certain amount of energy. The occupancy depends on the energy
of the state and the temperature and the chemical potential of the system.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. :math:`E_i > \mu`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics#Bose%E2%80%93Einstein_distribution>`__.
"""

from sympy import Eq, exp
from symplyphysics import (
    convert_to_float,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)

occupancy_of_state = Symbol("n_i", dimensionless)
"""
Occupancy of single-particle state :math:`i`.
"""

energy_of_state = clone_as_symbol(symbols.energy, subscript="i")
"""
:symbols:`energy` of single-particle state :math:`i`.
"""

total_chemical_potential = symbols.chemical_potential
"""
Total :symbols:`chemical_potential` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

law = Eq(
    occupancy_of_state, 1 / (exp(
    (energy_of_state - total_chemical_potential) / (quantities.boltzmann_constant * temperature)) - 1))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    energy_of_state_=energy_of_state,
    total_chemical_potential_=total_chemical_potential,
    temperature_=temperature,
)
@validate_output(occupancy_of_state)
def calculate_occupancy_of_state(
    energy_of_state_: Quantity,
    total_chemical_potential_: Quantity,
    temperature_: Quantity,
) -> float:
    result = law.rhs.subs({
        energy_of_state: energy_of_state_,
        total_chemical_potential: total_chemical_potential_,
        temperature: temperature_,
    })
    return convert_to_float(result)
