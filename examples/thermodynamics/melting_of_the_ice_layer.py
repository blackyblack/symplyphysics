#!/usr/bin/env python3

from sympy import solve, Symbol, Eq
from symplyphysics import print_expression, Quantity, prefixes, units, convert_to
from symplyphysics.core.symbols.celsius import to_kelvin_quantity, Celsius
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as energy_heating_law
from symplyphysics.laws.thermodynamics import energy_to_melt_from_mass as energy_melting_law
from symplyphysics.definitions import density_from_mass_volume as density_law
from symplyphysics.laws.thermodynamics import sum_of_heat_transfer_is_zero as thermodinamics_law_1

# The 4.2 cm thick ice layer has a temperature of 0 °C.
# What is the minimum thickness of the water layer at a temperature of 306 K
# that needs to be poured on ice so that it all melts?

layer_thickness_of_ice = Symbol("height_of_ice")
temperature_of_ice = Symbol("temperature_of_ice")
temperature_of_water = Symbol("temperature_of_water")

density_of_water = Symbol("density_of_water")
density_of_ice = Symbol("density_of_ice")
specific_heat_heating_of_water = Symbol("specific_heat_heating_of_water")
specific_heat_melting_of_ice = Symbol("specific_heat_melting_of_ice")

layer_thickness_of_water = Symbol("height_of_water")

# In this example, the surface area of ice is equal to the surface area of water.
# The thicknesses of the ice and water layers are related to the volumes by the formula:
# V = h * S
volume = Symbol("volume")
height = Symbol("height")
area = Symbol("area")
volume_equation = Eq(volume, height * area)

volume_of_ice_value = volume_equation.subs({
    height: layer_thickness_of_ice
}).rhs
volume_of_water_value = volume_equation.subs({
    height: layer_thickness_of_water
}).rhs

density_of_ice_equation = density_law.definition.subs({
    density_law.density: density_of_ice,
    density_law.volume: volume_of_ice_value
})
mass_of_ice = solve(density_of_ice_equation, density_law.mass, dict=True)[0][density_law.mass]

density_of_water_equation = density_law.definition.subs({
    density_law.density: density_of_water,
    density_law.volume: volume_of_water_value
})
mass_of_water = solve(density_of_water_equation, density_law.mass, dict=True)[0][density_law.mass]

energy_for_melting_ice_value = energy_melting_law.law.subs({
    energy_melting_law.mass_of_matter: mass_of_ice,
    energy_melting_law.specific_heat_melting: specific_heat_melting_of_ice
}).rhs
energy_from_cooling_water_value = energy_heating_law.law.subs({
    energy_heating_law.body_mass: mass_of_water,
    energy_heating_law.specific_heat_capacity: specific_heat_heating_of_water,
    energy_heating_law.temperature_origin: temperature_of_water,
    energy_heating_law.temperature_end: temperature_of_ice
}).rhs

heat_balance_equation = thermodinamics_law_1.law.subs({
    thermodinamics_law_1.amounts_energy: (energy_for_melting_ice_value, energy_from_cooling_water_value)
}).doit()

print(print_expression(heat_balance_equation))

height_of_water_solve = solve(heat_balance_equation, layer_thickness_of_water, dict=True)[0][layer_thickness_of_water]
answer = Eq(layer_thickness_of_water, height_of_water_solve)
print(f"Total equation:\n{print_expression(answer)}")

height_of_water_value = height_of_water_solve.subs({
    layer_thickness_of_ice: Quantity(4.2 * prefixes.centi * units.meters),
    temperature_of_ice: to_kelvin_quantity(Celsius(0)),
    temperature_of_water: Quantity(306 * units.kelvins) ,
    density_of_water: Quantity(1000 * units.kilograms / (units.meter ** 3)),
    density_of_ice: Quantity(900 * units.kilograms / (units.meter ** 3)),
    specific_heat_heating_of_water: Quantity(4200 * units.joules / (units.kilogram * units.kelvin)),
    specific_heat_melting_of_ice: Quantity(330 * prefixes.kilo * units.joules / units.kilogram)
})
print(f"Need height of water is: {convert_to(Quantity(height_of_water_value), units.meters)}")
