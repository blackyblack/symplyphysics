from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

from symplyphysics.laws.thermodynamics import pressure_from_temperature_and_volume as ideal_gas_law
from symplyphysics.laws.thermodynamics import average_kinetic_energy_of_molecules_from_temperature as kinetic_energy
from symplyphysics.laws.chemistry import avogadro_number_from_mole_count as avogadro_number
from symplyphysics.definitions import volume_number_density
# Description

## Ideal gas equation: p = 2/3 * n * E
## Where:
## p is pressure
## n is concentration of molecules
## E is the average kinetic energy of gas molecules

## Conditions
## The gas must be ideal

average_kinetic_energy = Symbol("average_kinetic_energy", units.energy)
molecules_concentration = Symbol("molecules_concentration", 1 / units.volume)
pressure = Symbol("pressure", units.pressure)

law = Eq(pressure, 2 / 3 * molecules_concentration * average_kinetic_energy)

## Proof

temperature_eq = kinetic_energy.law.subs(
    {kinetic_energy.average_kinetic_energy: average_kinetic_energy})

derived_temperature = solve(temperature_eq, kinetic_energy.temperature,
    dict=True)[0][kinetic_energy.temperature]

particles_number_eq = volume_number_density.definition.subs({
    volume_number_density.volume: ideal_gas_law.volume,
    volume_number_density.number_density: molecules_concentration
})

derived_particles_number = solve(particles_number_eq, volume_number_density.objects,
    dict=True)[0][volume_number_density.objects]

mole_count_eq = avogadro_number.law.subs(
    {avogadro_number.particles_count: derived_particles_number})

derived_mole_count = solve(mole_count_eq, avogadro_number.mole_count,
    dict=True)[0][avogadro_number.mole_count]

derived_pressure = ideal_gas_law.law.subs({
    ideal_gas_law.temperature: derived_temperature,
    ideal_gas_law.units.molar_gas_constant: units.boltzmann * units.avogadro,
    ideal_gas_law.mole_count: derived_mole_count
})

assert derived_pressure.rhs == law.rhs


def print_law() -> str:
    return print_expression(law)


@validate_input(molecules_concentration_=molecules_concentration,
    average_kinetic_energy_=average_kinetic_energy)
@validate_output(pressure)
def calculate_pressure(molecules_concentration_: Quantity,
    average_kinetic_energy_: Quantity) -> Quantity:
    result_expr = solve(law, pressure, dict=True)[0][pressure]
    result_pressure = result_expr.subs({
        molecules_concentration: molecules_concentration_,
        average_kinetic_energy: average_kinetic_energy_,
    })
    return Quantity(result_pressure)
