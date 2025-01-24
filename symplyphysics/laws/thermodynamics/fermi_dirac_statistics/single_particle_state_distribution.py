"""
Single-particle state distribution
==================================

For a system of identical fermions in thermodynamic equilibrium, the average number of fermions
in a single-particle state :math:`i` is given by the Fermi—Dirac distribution.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Notes:**

#. If the energy states are degenerate, i.e. two or more particles are on the same energy level,
   the average number of fermions can be found by multiplying by the degeneracy :math:`g_i` of
   the energy level.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Fermi%E2%80%93Dirac_statistics#Fermi%E2%80%93Dirac_distribution>`__.
"""

from sympy import Eq, exp
from symplyphysics import (
    dimensionless,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
    quantities,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

occupancy_of_state = SymbolNew("N_i", dimensionless)
"""
Occupancy of a single-particle state :math:`i`.
"""

energy_of_state = clone_as_symbol(symbols.energy, subscript="i")
"""
:symbols:`energy` of state :math:`i`.
"""

total_chemical_potential = symbols.chemical_potential
r"""
Total :symbols:`chemical_potential` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

law = Eq(
    occupancy_of_state, 1 / (exp((energy_of_state - total_chemical_potential) /
    (quantities.boltzmann_constant * temperature)) + 1))
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
