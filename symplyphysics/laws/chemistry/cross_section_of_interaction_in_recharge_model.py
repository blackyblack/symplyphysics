from sympy import Eq, nsolve, pi, log, sqrt
from sympy.physics.units import elementary_charge, boltzmann_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    clone_symbol,
    symbols
)

# Description
## The effective cross section is a physical quantity characterizing the probability of transition of a system of
## two interacting particles to a certain final state, a quantitative characteristic of the acts of collision of
## particles of a stream hitting a target with target particles. The effective cross-section has the dimension of the area.

## Law is: g = pi * a0^2 * (Uh / Ui) * (ln(sqrt(3 * k * T / m) * (g * P * m / (2 * k * T * q * E)) * sqrt(Ui / Uh)))^2, where
## g - the cross-sectional area of the interaction of particles,
## Ui - the ionization energy of atoms,
## T - temperature,
## m - mass of atom,
## P - pressure
## E - electric intensity,
## Uh - hydrogen ionization energy,
## a0 - bohr radius,
## k - boltzmann constant,
## q - elementary charge.

cross_sectional_area_of_interaction = Symbol("cross_sectional_area_of_interaction", units.area)

ionization_energy = Symbol("ionization_energy", units.energy)
mass_of_atom = clone_symbol(symbols.basic.mass, "mass_of_atom")
pressure = Symbol("pressure", units.pressure)
temperature = clone_symbol(symbols.thermodynamics.temperature, "temperature")
electric_intensity = Symbol("electric_intensity", units.voltage / units.length)

bohr_radius = Quantity(0.529e-10 * units.meter)
hydrogen_ionization_energy = Quantity(13.6 * units.electronvolt)

expression_1 = sqrt(ionization_energy / hydrogen_ionization_energy)
expression_2 = sqrt(3 * boltzmann_constant * temperature / mass_of_atom)
expression_3 = cross_sectional_area_of_interaction * pressure * mass_of_atom
expression_4 = 2 * boltzmann_constant * temperature * elementary_charge * electric_intensity
expression_5 = hydrogen_ionization_energy / ionization_energy

law = Eq(cross_sectional_area_of_interaction, pi * bohr_radius**2 * expression_5 * (log(expression_2 * (expression_3 / expression_4) * expression_1))**2)


@validate_input(ionization_energy_=ionization_energy,
    mass_of_atom_=mass_of_atom,
    pressure_=pressure,
    temperature_=temperature,
    electric_intensity_=electric_intensity)
@validate_output(cross_sectional_area_of_interaction)
def calculate_cross_sectional_area_of_interaction(ionization_energy_: Quantity,
    mass_of_atom_: Quantity, pressure_: Quantity,
    temperature_: Quantity, electric_intensity_: Quantity) -> Quantity:
    # nsolve() only works with numerical equations
    applied_law = law.subs({
        ionization_energy: ionization_energy_.scale_factor,
        mass_of_atom: mass_of_atom_.scale_factor,
        pressure: pressure_.scale_factor,
        temperature: temperature_.scale_factor,
        electric_intensity: electric_intensity_.scale_factor,
        bohr_radius: bohr_radius.scale_factor,
        hydrogen_ionization_energy: hydrogen_ionization_energy.scale_factor,
        elementary_charge: elementary_charge.scale_factor,
        boltzmann_constant: boltzmann_constant.scale_factor
    })
    result_expr = nsolve(applied_law, cross_sectional_area_of_interaction, 1)
    return Quantity(result_expr, dimension=units.area)
