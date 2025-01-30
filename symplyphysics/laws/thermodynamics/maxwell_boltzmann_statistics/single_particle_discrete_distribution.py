r"""
Discrete distribution
=====================

Maxwell—Boltzmann distribution can be written as a discrete distribution of a single particle's
discrete energy spectrum. Maxwell-Boltzmann statistics gives the average number of particles
found in a given single-particle microstate.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. Particles do not interact and are classical.
#. The system is in thermal equilibrium.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_statistics>`__.
"""

from sympy import Eq, exp, solve
from symplyphysics import (
    dimensionless,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
    quantities,
    clone_as_symbol,
)

occupancy_of_state = SymbolNew("N_i", dimensionless)
"""
Occupancy of, or expected number of particles in, the single-particle microstate :math:`i`.
"""

particle_count = symbols.particle_count
"""
Total :symbols:`particle_count` of the system.
"""

energy_of_state = clone_as_symbol(symbols.energy, subscript="i")
"""
:symbols:`energy` of single-particle microstate :math:`i`.
"""

equilibrium_temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

single_particle_partition_function = symbols.partition_function
"""
Single-particle :symbols:`partition_function`, which acts as a normalizing factor of the distribution.
"""

law = Eq(occupancy_of_state, (particle_count / single_particle_partition_function) *
    exp(-1 * energy_of_state / (quantities.boltzmann_constant * equilibrium_temperature)))
"""
:laws:symbol::

:laws:latex::
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
