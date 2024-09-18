r"""
Free energy differential
========================

The fundamental thermodynamic relations are fundamental equations which demonstate how important
thermodynamic quantities depend on variables that are measurable experimentally.

**Notation:**

#. :math:`d` denotes an exact, path-independent differential.

**Notes:**

#. Temperature, volume, and particle count are so called natural variables of free energy as a
   thermodynamic potential.
#. For a system with more than one type of particles, the last term can be represented as a sum over all
   types of particles, i.e. :math:`\sum_i \mu_i \, d N_i`.

**Conditions:**

#. The system is in thermal equilibrium with its surroundings.
#. The system is composed of only one type of particles, i.e. the system is a pure substance.
"""

from sympy import Eq, symbols as sympy_symbols, Function as SymFunction, Symbol as SymSymbol
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
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import (
    helmholtz_free_energy_via_internal_energy as free_energy_def,
    internal_energy_differential,
)

free_energy_change = Symbol("free_energy_change", units.energy)
"""
Infinitesimal change in free energy of the system.

Symbol:
    :code:`dF`
"""

entropy = Symbol("entropy", units.energy / units.temperature)
"""
Enthalpy of the system.

Symbol:
    :code:`S`
"""

temperature_change = clone_as_symbol(symbols.temperature, display_symbol="dT")
"""
Infinitesimal change in :symbols:`temperature` of the system.
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
    free_energy_change,
    -1 * entropy * temperature_change - pressure * volume_change +
    chemical_potential * particle_count_change,
)
r"""
:code:`dF = -1 * S * dT - p * dV + mu * dN`

Latex:
    .. math::
        dF = -S \, dT - p \, dV + \mu \, dN
"""

# Derive from the definition of Helmholtz free energy and internal energy differential

_internal_energy_sym = sympy_symbols("internal_energy", cls=SymFunction)
_entropy_change = SymSymbol("entropy_change")
_temperature = SymSymbol("temperature")
_volume = SymSymbol("volume")
_particle_count = SymSymbol("particle_count")

_internal_energy = _internal_energy_sym(entropy, _volume, _particle_count)

_free_energy = free_energy_def.law.rhs.subs({
    free_energy_def.internal_energy: _internal_energy,
    free_energy_def.temperature: _temperature,
    free_energy_def.entropy: entropy,
})

# The differential of free energy can be found by adding up the products of
# the partial derivative of free energy with respect to a thermodynamic quantity
# and the infinitesimal change of that quantity
_free_energy_change = (_free_energy.diff(_temperature) * temperature_change +
    _free_energy.diff(entropy) * _entropy_change + _free_energy.diff(_volume) * volume_change +
    _free_energy.diff(_particle_count) * particle_count_change)

_internal_energy_diff = internal_energy_differential.law.rhs.subs({
    internal_energy_differential.temperature: _temperature,
    internal_energy_differential.entropy_change: _entropy_change,
    internal_energy_differential.pressure: pressure,
    internal_energy_differential.volume_change: volume_change,
    internal_energy_differential.chemical_potential: chemical_potential,
    internal_energy_differential.particle_count_change: particle_count_change,
})

_internal_energy_diff_entropy = _internal_energy_diff.subs({
    _entropy_change: 1,
    volume_change: 0,
    particle_count_change: 0
})

_internal_energy_diff_volume = _internal_energy_diff.subs({
    _entropy_change: 0,
    volume_change: 1,
    particle_count_change: 0
})

_internal_energy_diff_particle_count = _internal_energy_diff.subs({
    _entropy_change: 0,
    volume_change: 0,
    particle_count_change: 1
})

_free_energy_change = _free_energy_change.subs({
    _internal_energy.diff(entropy): _internal_energy_diff_entropy,
    _internal_energy.diff(_volume): _internal_energy_diff_volume,
    _internal_energy.diff(_particle_count): _internal_energy_diff_particle_count,
})

assert expr_equals(_free_energy_change, law.rhs)


@validate_input(
    entropy_=entropy,
    temperature_change_=temperature_change,
    pressure_=pressure,
    volume_change_=volume_change,
    chemical_potential_=chemical_potential,
    particle_count_change_=particle_count_change,
)
@validate_output(free_energy_change)
# pylint: disable=too-many-arguments
def calculate_free_energy_change(
    entropy_: Quantity,
    temperature_change_: Quantity,
    pressure_: Quantity,
    volume_change_: Quantity,
    chemical_potential_: Quantity,
    particle_count_change_: int,
) -> Quantity:
    result = law.rhs.subs({
        entropy: entropy_,
        temperature_change: temperature_change_,
        pressure: pressure_,
        volume_change: volume_change_,
        chemical_potential: chemical_potential_,
        particle_count_change: particle_count_change_,
    })
    return Quantity(result)
