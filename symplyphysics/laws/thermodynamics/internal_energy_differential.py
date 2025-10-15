"""
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
   types of particles, i.e. :math:`\\sum_i \\mu_i \\, d N_i`.

**Conditions:**

#. The system is in thermal equilibrium with its surroundings
#. The system is composed of only one type of particles, i.e. the system is a pure substance.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Fundamental_thermodynamic_relation>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.operations.symbolic import ExactDifferential
from symplyphysics.laws.thermodynamics import (
    entropy_change_in_reversible_process as _second_law,
    infinitesimal_work_in_quasistatic_process as _work_law,
    internal_energy_change_via_heat_and_work as _first_law,
)

internal_energy_change = ExactDifferential(symbols.internal_energy)
"""
Infinitesimal change in :symbols:`internal_energy` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

entropy_change = ExactDifferential(symbols.entropy)
"""
Infinitesimal change in :symbols:`entropy` of the system.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` inside the system.
"""

volume_change = ExactDifferential(symbols.volume)
"""
Infinitesimal change in :symbols:`volume` of the system.
"""

chemical_potential = symbols.chemical_potential
"""
:symbols:`chemical_potential` of the system.
"""

particle_count_change = ExactDifferential(symbols.particle_count)
"""
Infinitesimal change in the number of particles in the system.
"""

law = Eq(
    internal_energy_change,
    temperature * entropy_change - pressure * volume_change +
    chemical_potential * particle_count_change,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from the first and second laws of thermodynamics
# TODO: refactor proof to take account of changes in the number of particles

_heat_supplied_to_system = solve(_second_law.law, _second_law.heat_supplied_to_system)[0].subs({
    _second_law.entropy_change: entropy_change,
    _second_law.common_temperature: temperature,
})

_work_done_by_system = _work_law.law.rhs.subs({
    _work_law.pressure: pressure,
    _work_law.volume_change: volume_change,
})

_internal_energy_change = _first_law.law.rhs.subs({
    _first_law.heat_supplied_to_system: _heat_supplied_to_system,
    _first_law.work_done_by_system: _work_done_by_system,
})

_expr_from_law = law.rhs.subs({particle_count_change: 0})

assert expr_equals(_expr_from_law, _internal_energy_change)


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


# UNIQUE_LAW_ID: 60
