r"""
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
"""

from sympy import Eq, exp
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    quantities,
    convert_to_float,
)

occupancy_of_state = Symbol("occupancy_of_state", dimensionless)
r"""
Occupancy of a single-particle state :math:`i`.

Symbol:
    :code:`N_i`

Latex:
    :math:`N_i`
"""

energy_of_state = Symbol("energy_of_state", units.energy)
r"""
Energy of state :math:`i`.

Symbol:
    :code:`E_i`

Latex:
    :math:`E_i`
"""

total_chemical_potential = Symbol("total_chemical_potential", units.energy)
r"""
Total chemical potential of the system.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

law = Eq(
    occupancy_of_state, 1 / (exp(
    (energy_of_state - total_chemical_potential) / (quantities.boltzmann_constant * temperature)) + 1))
r"""
:code:`N_i = 1 / (exp((E_i - mu) / (k_B * T)) + 1)`

Latex:
    .. math::
        N_i = \frac{1}{\exp{\left( \frac{E_i - \mu}{k_\text{B} T} \right)} + 1}
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
