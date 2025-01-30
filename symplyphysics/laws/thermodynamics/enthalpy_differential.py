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

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Fundamental_thermodynamic_relation>`__.
"""

from sympy import Eq, Function as SymFunction, Symbol as SymSymbol, symbols as sympy_symbols
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import (
    enthalpy_is_internal_energy_plus_pressure_energy as enthalpy_def,
    internal_energy_differential,
)

enthalpy_change = clone_as_symbol(symbols.enthalpy, display_symbol="dH", display_latex="dH")
"""
Infinitesimal change in :symbols:`enthalpy` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

entropy_change = clone_as_symbol(symbols.entropy, display_symbol="dS", display_latex="dS")
"""
Infinitesimal change in :symbols:`entropy` of the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

pressure_change = clone_as_symbol(symbols.pressure, display_symbol="dp", display_latex="dp")
"""
Infinitesimal change in :symbols:`pressure` inside the system.
"""

chemical_potential = symbols.chemical_potential
"""
:symbols:`chemical_potential` of the system.
"""

particle_count_change = clone_as_symbol(symbols.particle_count, display_symbol="dN", display_latex="dN")
"""
Infinitesimal change in the :symbols:`particle_count` of the system.
"""

law = Eq(
    enthalpy_change,
    temperature * entropy_change + volume * pressure_change +
    chemical_potential * particle_count_change,
)
"""
:laws:symbol::

:laws:latex::
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
def calculate_enthalpy_change(
    temperature_: Quantity,
    entropy_change_: Quantity,
    volume_: Quantity,
    pressure_change_: Quantity,
    chemical_potential_: Quantity,
    particle_count_change_: int,
) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    result = law.rhs.subs({
        temperature: temperature_,
        entropy_change: entropy_change_,
        volume: volume_,
        pressure_change: pressure_change_,
        chemical_potential: chemical_potential_,
        particle_count_change: particle_count_change_,
    })
    return Quantity(result)
