r"""
Single particle state distribution
==================================

Occupancy of a single-particle state of bosons is the probability of a single boson
to occupy a state with a certain amount of energy. The occupancy depends on the energy
of the state and the temperature and the chemical potential of the system.

**Conditions:**

#. :math:`E_i > \mu`.
"""

from sympy import Eq, exp
from symplyphysics import (
    convert_to_float,
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

occupancy_of_state = Symbol("occupancy_of_state", dimensionless)
r"""
Occupancy of single-particle state :math:`i`.

Symbol:
    :code:`n_i`

Latex:
    :math:`n_i`
"""

energy_of_state = Symbol("energy_of_state", units.energy)
r"""
Energy of single-particle state :math:`i`.

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
    (energy_of_state - total_chemical_potential) / (units.boltzmann * temperature)) - 1))
r"""
:code:`n_i = 1 / (exp((E_i - mu) / (k_B * T)) - 1)`

Latex:
    .. math::
        n_i = \frac{1}{\exp{\frac{E_i - \mu}{k_\text{B} T}} - 1}
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
