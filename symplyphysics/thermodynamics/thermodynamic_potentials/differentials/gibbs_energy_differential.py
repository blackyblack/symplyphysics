"""
Gibbs energy differential
=========================

The fundamental thermodynamic relations are fundamental equations which demonstate how important
thermodynamic quantities depend on variables that are measurable experimentally.

**Notation:**

#. :math:`d` denotes an exact, path-independent differential.

**Notes:**

#. Temperature, pressure, and particle count are so called natural variables of Gibbs energy as a
   thermodynamic potential.
#. For a system with more than one type of particles, the last term can be represented as a sum over all
   types of particles, i.e. :math:`\\sum_i \\mu_i \\, d N_i`.

**Conditions:**

#. The system is in thermal equilibrium with its surroundings.
#. The system is composed of only one type of particles, i.e. the system is a pure substance.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Fundamental_thermodynamic_relation>`__.
"""

from sympy import Eq, symbols as sympy_symbols, Function as SymFunction, Symbol as SymSymbol
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.thermodynamics.thermodynamic_potentials.definitions import gibbs_energy_via_enthalpy as gibbs_energy_def
from symplyphysics.thermodynamics.thermodynamic_potentials.differentials import enthalpy_differential

gibbs_energy_change = clone_as_symbol(symbols.gibbs_energy, display_symbol="dG", display_latex="dG")
"""
Infinitesimal change in :symbols:`gibbs_energy` of the system.
"""

entropy = symbols.entropy
"""
:symbols:`entropy` of the system.
"""

temperature_change = clone_as_symbol(symbols.temperature, display_symbol="dT", display_latex="dT")
"""
Infinitesimal change in :symbols:`temperature` of the system.
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

particle_count_change = clone_as_symbol(symbols.particle_count,
    display_symbol="dN",
    display_latex="dN")
"""
Infinitesimal change in the :symbols:`particle_count` of the system.
"""

law = Eq(
    gibbs_energy_change,
    -1 * entropy * temperature_change + volume * pressure_change +
    chemical_potential * particle_count_change,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from the definition of Gibbs energy and enthalpy differential

_enthalpy_sym = sympy_symbols("enthalpy", cls=SymFunction)
_entropy_change = SymSymbol("entropy_change")
_temperature = SymSymbol("temperature")
_pressure = SymSymbol("pressure")
_particle_count = SymSymbol("particle_count")

_enthalpy = _enthalpy_sym(entropy, _pressure, _particle_count)

_gibbs_energy = gibbs_energy_def.law.rhs.subs({
    gibbs_energy_def.enthalpy: _enthalpy,
    gibbs_energy_def.temperature: _temperature,
    gibbs_energy_def.entropy: entropy,
})

# The differential of Gibbs energy can be found by adding up the products of
# the partial derivative of Gibbs energy with respect to a thermodynamic quantity
# and the infinitesimal change of that quantity
_gibbs_energy_change = (_gibbs_energy.diff(_temperature) * temperature_change +
    _gibbs_energy.diff(entropy) * _entropy_change +
    _gibbs_energy.diff(_pressure) * pressure_change +
    _gibbs_energy.diff(_particle_count) * particle_count_change)

_enthalpy_diff = enthalpy_differential.law.rhs.subs({
    enthalpy_differential.temperature: _temperature,
    enthalpy_differential.entropy_change: _entropy_change,
    enthalpy_differential.volume: volume,
    enthalpy_differential.pressure_change: pressure_change,
    enthalpy_differential.chemical_potential: chemical_potential,
    enthalpy_differential.particle_count_change: particle_count_change,
})

_enthalpy_diff_entropy = _enthalpy_diff.subs({
    _entropy_change: 1,
    pressure_change: 0,
    particle_count_change: 0
})

_enthalpy_diff_pressure = _enthalpy_diff.subs({
    _entropy_change: 0,
    pressure_change: 1,
    particle_count_change: 0
})

_enthalpy_diff_particle_count = _enthalpy_diff.subs({
    _entropy_change: 0,
    pressure_change: 0,
    particle_count_change: 1
})

_gibbs_energy_change = _gibbs_energy_change.subs({
    _enthalpy.diff(entropy): _enthalpy_diff_entropy,
    _enthalpy.diff(_pressure): _enthalpy_diff_pressure,
    _enthalpy.diff(_particle_count): _enthalpy_diff_particle_count,
})

assert expr_equals(_gibbs_energy_change, law.rhs)


@validate_input(
    entropy_=entropy,
    temperature_change_=temperature_change,
    volume_=volume,
    pressure_change_=pressure_change,
    chemical_potential_=chemical_potential,
    particle_count_change_=particle_count_change,
)
@validate_output(gibbs_energy_change)
def calculate_gibbs_energy_change(
    entropy_: Quantity,
    temperature_change_: Quantity,
    volume_: Quantity,
    pressure_change_: Quantity,
    chemical_potential_: Quantity,
    particle_count_change_: int,
) -> Quantity:
    result = law.rhs.subs({
        entropy: entropy_,
        temperature_change: temperature_change_,
        volume: volume_,
        pressure_change: pressure_change_,
        chemical_potential: chemical_potential_,
        particle_count_change: particle_count_change_,
    })
    return Quantity(result)
