r"""
Discrete distribution
=====================

Maxwell—Boltzmann distribution can be written as a discrete distribution of a single particle's
discrete energy spectrum. Maxwell-Boltzmann statistics gives the average number of particles
found in a given single-particle microstate.

**Notation:**

#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.

**Conditions:**

#. Particles do not interact and are classical.
#. The system is in thermal equilibrium.
"""

from sympy import Eq, exp, solve
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

occupancy_of_state = Symbol("occupancy_of_state", dimensionless)
r"""
Occupancy of, or expected number of particles in, the single-particle microstate
:math:`i`.

Symbol:
    :code:`N_i`

Latex:
    :math:`N_i`
"""

particle_count = Symbol("particle_count", dimensionless)
"""
Total number of particles in the system.

Symbol:
    :code:`N`
"""

energy_of_state = Symbol("energy_of_state", units.energy)
r"""
Energy of single-particle microstate :math:`i`.

Symbol:
    :code:`E_i`

Latex:
    :math:`E_i`
"""

equilibrium_temperature = symbols.temperature
"""
:attr:`~symplyphysics.symbols.temperature` of the system.
"""

single_particle_partition_function = Symbol("single_particle_partition_function", dimensionless)
"""
Single-particle partition function, which acts as a normalizing factor of the distribution.

Symbol:
    :code:`Z`
"""

law = Eq(occupancy_of_state, (particle_count / single_particle_partition_function) *
    exp(-1 * energy_of_state / (units.boltzmann_constant * equilibrium_temperature)))
r"""
:code:`N_i = (N / Z) * exp(-1 * E_i / (k_B * T))`

Latex:
    .. math::
        N_i = \frac{N}{Z} \exp \left( - \frac{E_i}{k_\text{B} T} \right)
"""


@validate_input(
    total_particle_count_=particle_count,
    energy_of_microstate_=energy_of_state,
    equilibrium_temperature_=equilibrium_temperature,
    single_particle_partition_function_=single_particle_partition_function,
)
@validate_output(occupancy_of_state)
def calculate_particle_count_in_microstate(
    total_particle_count_: int,
    energy_of_microstate_: Quantity,
    equilibrium_temperature_: Quantity,
    single_particle_partition_function_: float,
) -> Quantity:
    expr = solve(law, occupancy_of_state)[0]

    result = expr.subs({
        particle_count: total_particle_count_,
        energy_of_state: energy_of_microstate_,
        equilibrium_temperature: equilibrium_temperature_,
        single_particle_partition_function: single_particle_partition_function_,
    })
    return Quantity(result)
