"""
Cross section of interaction in recharge model
==============================================

The effective cross section is a physical quantity characterizing the probability of
transition of a system of two interacting particles to a certain final state, a
quantitative characteristic of the acts of collision of particles of a stream hitting a
target with target particles. The effective cross-section has the dimension of the area.

**Notation:**

#. :quantity_notation:`bohr_radius`.
#. :quantity_notation:`hydrogen_ionization_energy`.
#. :quantity_notation:`boltzmann_constant`.
#. :quantity_notation:`elementary_charge`.

..
    TODO: find link
    TODO: move to `ionization` folder?
"""

from sympy import Eq, nsolve, pi, log, sqrt, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.quantities import (
    bohr_radius,
    hydrogen_ionization_energy,
    boltzmann_constant,
    elementary_charge,
)
from symplyphysics.core.convert import evaluate_expression

cross_section = symbols.cross_section
"""
:symbols:`cross_section` of interaction of particles.
"""

ionization_energy = clone_as_symbol(symbols.energy,
    display_symbol="E_i",
    display_latex="E_\\text{i}")
"""
Ionization :symbols:`energy` of the particles.
"""

molecular_mass = symbols.mass
"""
:symbols:`mass` of a single gas particle.
"""

pressure = symbols.pressure
"""
Gas :symbols:`pressure`.
"""

temperature = symbols.temperature
"""
Gas :symbols:`temperature`.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

law = Eq(cross_section,
    (pi * bohr_radius**2) * (hydrogen_ionization_energy / ionization_energy) * log(
    sqrt(3 * boltzmann_constant * temperature / molecular_mass) *
    sqrt(ionization_energy / hydrogen_ionization_energy) *
    ((cross_section * pressure * molecular_mass) /
    (2 * boltzmann_constant * temperature * elementary_charge * electric_field_strength)))**2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(ionization_energy_=ionization_energy,
    mass_of_atom_=molecular_mass,
    pressure_=pressure,
    temperature_=temperature,
    electric_intensity_=electric_field_strength)
@validate_output(cross_section)
def calculate_cross_sectional_area_of_interaction(ionization_energy_: Quantity,
    mass_of_atom_: Quantity, pressure_: Quantity, temperature_: Quantity,
    electric_intensity_: Quantity) -> Quantity:
    # NOTE `nsolve` doesn't recognize `SymbolNew` instances, works fine with our old `Symbol` class though.
    cross_section_sym = SymSymbol("sigma")
    eqn = law.subs({
        cross_section: cross_section_sym,
        ionization_energy: ionization_energy_,
        molecular_mass: mass_of_atom_,
        pressure: pressure_,
        temperature: temperature_,
        electric_field_strength: electric_intensity_,
    })
    # nsolve() only works with numerical equations
    eqn = evaluate_expression(eqn)
    result_expr = nsolve(eqn, cross_section_sym, 1)
    return Quantity(result_expr, dimension=units.area)
