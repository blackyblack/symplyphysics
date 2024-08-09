r"""
Enthalpy differential
=====================

The fundamental thermodynamic relations are fundamental equations which demonstate how the important
thermodynamic quantities depend on variables that are measurable experimentally.

**Notation:**

#. :math:`d` denotes an exact, path-independent differential.

**Notes:**

#. Entropy, pressure, and particle count are so called natural variables of enthalpy as a
   thermodynamic potential.
#. For a system with more than one type of particles, the last term can be represented as a sum over all
   types of particles, i.e. :math:`\sum_i \mu_i \, d N_i`.

**Conditions:**

#. The system is in thermal equilibrium with its surroundings.
#. The system is composed of only one type of particles, i.e. the system is a pure substance.
"""

from sympy import Eq, Function as SymFunction, Symbol as SymSymbol, symbols as sympy_symbols
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
    enthalpy_is_internal_energy_plus_pressure_energy as enthalpy_def,
    internal_energy_differential,
)

enthalpy_change = Symbol("enthalpy_change", units.energy)
"""
Infinitesimal change in enthalpy of the system.

Symbol:
    :code:`dH`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the system.

Symbol:
    :code:`T`
"""

entropy_change = Symbol("entropy_change", units.energy / units.temperature)
"""
Infinitesimal change in entropy of the system.

Symbol:
    :code:`dS`
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
    enthalpy_change,
    temperature * entropy_change + volume * pressure_change +
    chemical_potential * particle_count_change,
)
r"""
:code:`dH = T * dS + V * dp + mu * dN`

Latex:
    .. math::
        dH = T \, dS + V \, dp + \mu \, dN
"""

# Derive from definition of enthalpy and internal energy differential

_internal_energy_sym = sympy_symbols("internal_energy", cls=SymFunction)
_temperature_change = SymSymbol("temperature_change")
_entropy = SymSymbol("entropy")
_pressure = SymSymbol("pressure")
_volume_change = SymSymbol("volume_change")
_particle_count = SymSymbol("particle_count")

_internal_energy = _internal_energy_sym(_entropy, volume, _particle_count)

_enthalpy = enthalpy_def.law.rhs.subs({
    enthalpy_def.internal_energy: _internal_energy,
    enthalpy_def.pressure: _pressure,
    enthalpy_def.volume: volume,
})

# The differential of enthalpy can be found by adding up the products of
# the partial derivative of enthalpy with respect to a thermodynamic quantity
# and the infinitesimal change of that quantity
_enthalpy_change = (_enthalpy.diff(temperature) * _temperature_change +
    _enthalpy.diff(_entropy) * entropy_change + _enthalpy.diff(_pressure) * pressure_change +
    _enthalpy.diff(volume) * _volume_change +
    _enthalpy.diff(_particle_count) * particle_count_change)

_internal_energy_diff = internal_energy_differential.law.rhs.subs({
    internal_energy_differential.temperature: temperature,
    internal_energy_differential.entropy_change: entropy_change,
    internal_energy_differential.pressure: _pressure,
    internal_energy_differential.volume_change: _volume_change,
    internal_energy_differential.chemical_potential: chemical_potential,
    internal_energy_differential.particle_count_change: particle_count_change,
})

_internal_energy_diff_entropy = _internal_energy_diff.subs({
    entropy_change: 1,
    _volume_change: 0,
    particle_count_change: 0
})

_internal_energy_diff_volume = _internal_energy_diff.subs({
    entropy_change: 0,
    _volume_change: 1,
    particle_count_change: 0
})

_internal_energy_diff_particle_count = _internal_energy_diff.subs({
    entropy_change: 0,
    _volume_change: 0,
    particle_count_change: 1
})

_enthalpy_change = _enthalpy_change.subs({
    _internal_energy.diff(_entropy): _internal_energy_diff_entropy,
    _internal_energy.diff(volume): _internal_energy_diff_volume,
    _internal_energy.diff(_particle_count): _internal_energy_diff_particle_count,
})

assert expr_equals(_enthalpy_change, law.rhs)


@validate_input(
    temperature_=temperature,
    entropy_change_=entropy_change,
    volume_=volume,
    pressure_change_=pressure_change,
    chemical_potential_=chemical_potential,
    particle_count_change_=particle_count_change,
)
@validate_output(enthalpy_change)
# pylint: disable=too-many-arguments
def calculate_enthalpy_change(
    temperature_: Quantity,
    entropy_change_: Quantity,
    volume_: Quantity,
    pressure_change_: Quantity,
    chemical_potential_: Quantity,
    particle_count_change_: int,
) -> Quantity:
    result = law.rhs.subs({
        temperature: temperature_,
        entropy_change: entropy_change_,
        volume: volume_,
        pressure_change: pressure_change_,
        chemical_potential: chemical_potential_,
        particle_count_change: particle_count_change_,
    })
    return Quantity(result)
