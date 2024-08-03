from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation as ideal_gas_law
from symplyphysics.laws.thermodynamics import average_kinetic_energy_of_ideal_gas_from_temperature as kinetic_energy
from symplyphysics.laws.chemistry import avogadro_number_from_mole_count as avogadro_number
from symplyphysics.definitions import number_density_is_number_of_objects_per_unit_volume
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

derived_temperature = solve(temperature_eq, kinetic_energy.equilibrium_temperature,
    dict=True)[0][kinetic_energy.equilibrium_temperature]

particles_number_eq = number_density_is_number_of_objects_per_unit_volume.definition.subs({
    number_density_is_number_of_objects_per_unit_volume.volume: ideal_gas_law.volume,
    number_density_is_number_of_objects_per_unit_volume.number_density: molecules_concentration
})

derived_particles_number = solve(particles_number_eq,
    number_density_is_number_of_objects_per_unit_volume.number_of_objects,
    dict=True)[0][number_density_is_number_of_objects_per_unit_volume.number_of_objects]

mole_count_eq = avogadro_number.law.subs(
    {avogadro_number.particles_count: derived_particles_number})

derived_mole_count = solve(mole_count_eq, avogadro_number.mole_count,
    dict=True)[0][avogadro_number.mole_count]

derived_pressure = ideal_gas_law.law.subs({
    ideal_gas_law.temperature: derived_temperature,
    ideal_gas_law.units.molar_gas_constant: units.boltzmann * units.avogadro,
    ideal_gas_law.mole_count: derived_mole_count
})

assert expr_equals(derived_pressure.rhs, law.rhs)


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
