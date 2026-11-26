"""
Interaction cross section in recharge model
===========================================

This law gives the :ref:`effective cross section <effective_cross_section>` for resonant
charge exchange, the process in which a singly charged ion captures an electron from a
neutral atom of the same species. The newly formed ion accelerates from rest in a uniform
electric field over one mean free path.
The model uses the speed reached at the end of this path and the ideal-gas neutral
number density :math:`p/(k_\\mathrm{B}T)`.

**Notation:**

#. :quantity_notation:`hydrogen_ionization_energy`.
#. :quantity_notation:`boltzmann_constant`.
#. :quantity_notation:`elementary_charge`.

**Links:**

#. V. T. Barchenko et al., `Modeling Processes in Plasma Electronics
   <https://techeiscatel.ru/media/books/LabPlasma.pdf#page=7>`__ (LETI, 2008),
   p. 7, Eq. (1.4).
"""

from sympy import E, Eq, Symbol, log, nsolve, sqrt
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    convert_to,
)
from symplyphysics.quantities import (
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

reference_area = Quantity(0.88e-20 * units.meter**2, display_symbol="sigma_0",
    display_latex="\\sigma_0")
"""Area scale for charge exchange, normalized to hydrogen."""

reference_speed = Quantity(2.2e6 * units.meter / units.second, display_symbol="v_0")
"""Speed scale for hydrogen, rescaled by the atom's ionization energy for other species."""

_energy_ratio = ionization_energy / hydrogen_ionization_energy
_area_scale = reference_area / _energy_ratio
_neutral_density = pressure / (boltzmann_constant * temperature)
_collision_speed = sqrt((2 * elementary_charge * electric_field_strength /
    (molecular_mass * _neutral_density * cross_section)))
_speed_ratio = reference_speed * sqrt(_energy_ratio) / _collision_speed

law = Eq(cross_section, _area_scale * log(_speed_ratio)**2)
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
    ## Demonstrate the low-speed solution when the law admits multiple roots.
    ## The largest cross section gives the shortest acceleration path and the 
    ## lowest collision speed. If this branch is unavailable, return the first available root.
    substitutions = {
        ionization_energy: ionization_energy_,
        molecular_mass: mass_of_atom_,
        pressure: pressure_,
        temperature: temperature_,
        electric_field_strength: electric_intensity_,
    }
    if any(value.scale_factor <= 0 for value in substitutions.values()):
        raise ValueError("Charge-exchange model inputs must be positive.")

    area_scale = Quantity(_area_scale.subs(substitutions))
    # Evaluate the speed ratio when the cross section equals area_scale.
    speed_ratio_at_scale = convert_to(
        _speed_ratio.subs(cross_section, area_scale).subs(substitutions), 1)

    # On the low-speed branch, log(speed_ratio_at_scale) minimum is 1 at
    # cross_section = area_scale, so speed_ratio_at_scale must be at least E.
    if speed_ratio_at_scale == E:
        return area_scale
    if speed_ratio_at_scale > E:
        # The logarithm of a positive number is at most half that number, so
        # sqrt(cross_section / area_scale) <= 2*log(speed_ratio_at_scale).
        initial_area = area_scale * (2 * log(speed_ratio_at_scale))**2
    else:
        # Start with a speed ratio <= 1/e and an area <= area_scale. Here the law's
        # right side is >= area_scale and decreases toward zero.
        initial_area = area_scale * min(1, 1 / (E * speed_ratio_at_scale))**2

    sigma = Symbol("sigma", positive=True)
    equation = evaluate_expression(law.subs(substitutions).subs(cross_section, sigma))
    result = nsolve(equation, sigma, evaluate_expression(initial_area).evalf(), solver="newton",
        prec=50)
    # Restore square metres after solving in numerical SI units.
    return Quantity(result * units.meter**2)
