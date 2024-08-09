r"""
Internal energy differential
============================

The *fundamental thermodynamic relations* are fundamental equations which demonstrate how important
thermodynamic quantities depend on variables that are measurable experimentally.

**Notation:**

#. :math:`d` denotes an exact, path-independent differential.

**Notes:**

#. Entropy, volume, and particle count are so called natural variables of internal energy as a
   thermodynamic potential.
#. For a system with more than one type of particles, the last term can be represented as a sum over all
   types of particles, i.e. :math:`\sum_i \mu_i \, d N_i`.

**Conditions:**

#. The system is in thermal equilibrium with its surroundings
#. The system is composed of only one type of particles, i.e. the system is a pure substance.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import (
    entropy_change_in_reversible_process as second_law,
    infinitesimal_work_in_quasistatic_process as work_law,
    internal_energy_change_via_heat_and_work as first_law,
)

internal_energy_change = Symbol("internal_energy_change", units.energy)
"""
Infinitesimal change in internal energy of the system.

Symbol:
    :code:`dU`
"""

temperature = symbols.thermodynamics.temperature
"""
Temperature of the system.

Symbol:
    :code:`T`
"""

entropy_change = Symbol("entropy_change", units.energy / units.temperature)
"""
Infinitesimal change in entropy of the system.

Symbol:
    :code:`dS`
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure inside the system.

Symbol:
    :code:`p`
"""

volume_change = Symbol("volume_change", units.volume)
"""
Infinitesimal change in volume of the system.

Symbol:
    :code:`dV`
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
    internal_energy_change,
    temperature * entropy_change - pressure * volume_change +
    chemical_potential * particle_count_change,
)
r"""
:code:`dU = T * dS - p * dV + mu * dN`

Latex:
    .. math::
        dU = T \, dS - p \, dV + \mu \, dN
"""

# Derive from the first and second laws of thermodynamics
# TODO: refactor proof to take account of changes in the number of particles

_heat_supplied_to_system = solve(second_law.law,
    second_law.heat)[0].subs({
    second_law.entropy_change: entropy_change,
    second_law.common_temperature: temperature,
    })

_work_done_by_system = work_law.law.rhs.subs({
    work_law.pressure: pressure,
    work_law.infinitesimal_volume_change: volume_change,
})

_internal_energy_change = first_law.law.rhs.subs({
    first_law.heat_supplied_to_system: _heat_supplied_to_system,
    first_law.work_done_by_system: _work_done_by_system,
})

_expr_from_law = law.rhs.subs({particle_count_change: 0})

assert expr_equals(_expr_from_law, _internal_energy_change)


#pylint: disable=too-many-arguments
@validate_input(
    temperature_=temperature,
    entropy_change_=entropy_change,
    pressure_=pressure,
    volume_change_=volume_change,
    chemical_potential_=chemical_potential,
    particle_count_change_=particle_count_change,
)
@validate_output(internal_energy_change)
def calculate_internal_energy_change(
    temperature_: Quantity,
    entropy_change_: Quantity,
    pressure_: Quantity,
    volume_change_: Quantity,
    chemical_potential_: Quantity,
    particle_count_change_: int,
) -> Quantity:
    result = law.rhs.subs({
        temperature: temperature_,
        entropy_change: entropy_change_,
        pressure: pressure_,
        volume_change: volume_change_,
        chemical_potential: chemical_potential_,
        particle_count_change: particle_count_change_,
    })
    return Quantity(result)
