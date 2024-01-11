#!/usr/bin/env python3

from sympy import solve, Symbol, Eq
from symplyphysics import print_expression, Quantity, prefixes, units, convert_to
from symplyphysics.core.symbols.celsius import to_kelvin_quantity, Celsius
from symplyphysics.core.symbols.symbols import tuple_of_symbols
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as energy_heating_law
from symplyphysics.laws.thermodynamics import energy_to_melt_from_mass as energy_melting_law
from symplyphysics.definitions import density_from_mass_volume as density_law
from symplyphysics.laws.thermodynamics import sum_of_heat_transfer_is_zero as thermodinamics_law_1

# Example from https://easyfizika.ru/zadachi/termodinamika/vannu-emkostyu-100-litrov-neobhodimo-zapolnit-vodoj-imeyushhej-temperaturu-30-c/
# A bath with a capacity of 100 liters must be filled
# with water having a temperature of 30 °C,
# using water with a temperature of 80 °C and ice with a temperature of -20 °C.
# Find the mass of ice that needs to be put in the bath.

volume_of_bath = Symbol("volume_of_bath")
temperature_of_hot_water = Symbol("temperature_of_hot_water")
temperature_of_ice = Symbol("temperature_of_ice")
temperature_end = Symbol("temperature_end")

density_of_water = Symbol("density_of_water")
specific_heat_heating_water = Symbol("specific_heat_heating_water")
specific_heat_heating_ice = Symbol("specific_heat_heating_ice")
specific_heat_melting_ice = Symbol("specific_heat_melting_ice")
temperature_melt_ice = Symbol("temperature_melt_ice")

mass_of_ice = Symbol("mass_of_ice")
mass_of_hot_water = Symbol("mass_of_hot_water")

density_of_water_equation = density_law.definition.subs({
    density_law.volume: volume_of_bath,
    density_law.density: density_of_water
})
mass_of_all_water = solve(density_of_water_equation, density_law.mass,
    dict=True)[0][density_law.mass]

# the mass of all the water filling the bath consists of the mass of hot water
# that was in the bathroom initially, and the mass of water of melted ice
# mass_all_water = mass_of_hot_water + mass_of_ice
mass_of_all_water_equation = Eq(mass_of_all_water, mass_of_ice + mass_of_hot_water)
mass_of_hot_water_value = solve(mass_of_all_water_equation, mass_of_hot_water,
    dict=True)[0][mass_of_hot_water]

energy_cooling_hot_water = energy_heating_law.law.subs({
    energy_heating_law.specific_heat_capacity: specific_heat_heating_water,
    energy_heating_law.body_mass: mass_of_hot_water_value,
    energy_heating_law.temperature_origin: temperature_of_hot_water,
    energy_heating_law.temperature_end: temperature_end
})

energy_to_heating_ice_equation = energy_heating_law.law.subs({
    energy_heating_law.specific_heat_capacity: specific_heat_heating_ice,
    energy_heating_law.body_mass: mass_of_ice,
    energy_heating_law.temperature_origin: temperature_of_ice,
    energy_heating_law.temperature_end: temperature_melt_ice
})

energy_to_melt_ice_equation = energy_melting_law.law.subs({
    energy_melting_law.specific_heat_melting: specific_heat_melting_ice,
    energy_melting_law.mass_of_matter: mass_of_ice
})

energy_to_heat_melted_ice_equation = energy_heating_law.law.subs({
    energy_heating_law.specific_heat_capacity: specific_heat_heating_water,
    energy_heating_law.body_mass: mass_of_ice,
    energy_heating_law.temperature_origin: temperature_melt_ice,
    energy_heating_law.temperature_end: temperature_end
})

thermodinamics_law_1_equation = thermodinamics_law_1.law.subs({
    thermodinamics_law_1.amounts_energy:
    (energy_cooling_hot_water.rhs, energy_to_heating_ice_equation.rhs,
    energy_to_melt_ice_equation.rhs, energy_to_heat_melted_ice_equation.rhs)
}).doit()

mass_of_ice_equation = solve(thermodinamics_law_1_equation, mass_of_ice, dict=True)[0][mass_of_ice]
answer = Eq(mass_of_ice, mass_of_ice_equation)
print(f"Total equation:\n{print_expression(answer)}")

mass_of_ice_value_kg = mass_of_ice_equation.subs({
    volume_of_bath: Quantity(100 * units.liters),
    temperature_of_hot_water: to_kelvin_quantity(Celsius(80)),
    temperature_of_ice: to_kelvin_quantity(Celsius(-20)),
    temperature_end: to_kelvin_quantity(Celsius(30)),
    density_of_water: Quantity(1000 * units.kilograms / (units.meter**3)),
    specific_heat_heating_water: Quantity(4200 * units.joules / (units.kilogram * units.kelvin)),
    specific_heat_heating_ice: Quantity(2100 * units.joules / (units.kilogram * units.kelvin)),
    specific_heat_melting_ice: Quantity(330 * prefixes.kilo * units.joules / units.kilogram),
    temperature_melt_ice: to_kelvin_quantity(Celsius(0))
})
print(f"Mass of ice is: {convert_to(Quantity(mass_of_ice_value_kg), units.kilograms)} kg")
