#!/usr/bin/env python3

from sympy import solve, Symbol, Eq
from symplyphysics.core.symbols.celsius import to_kelvin_quantity, Celsius
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as energy_heating_law
from symplyphysics import print_expression, Quantity, prefixes, units, convert_to
from symplyphysics.laws.electricity import power_factor_from_active_and_full_power as efficiency_law
from symplyphysics.laws.thermodynamics import energy_from_combustion as combustion_energy_law

# Example from https://easyfizika.ru/zadachi/termodinamika/chtoby-nagret-1-8-kg-vody-ot-18-c-do-kipeniya-na-gorelke-s-kpd-25-potrebovalos/
# It took 92 g of fuel to heat 1.8 kg of water from 18 ° C to boiling on a burner
# with an efficiency of 25%. Find the specific heat of combustion of the fuel.

efficiency_factor = Symbol("efficiency_factor")
mass_of_fuel = Symbol("mass_of_fuel")
mass_of_liquid = Symbol("mass_of_liquid")
temperature_start = Symbol("temperature_start")
temperature_end = Symbol("temperature_end")

specific_heat_heating = Symbol("specific_heat_heating")

specific_heat_combustion = Symbol("specific_heat_of_combustion")

energy_heating_value = energy_heating_law.law.subs({
    energy_heating_law.specific_heat_capacity: specific_heat_heating,
    energy_heating_law.temperature_origin: temperature_start,
    energy_heating_law.temperature_end: temperature_end,
    energy_heating_law.body_mass: mass_of_liquid
}).rhs

energy_combustion_value = combustion_energy_law.law.subs({
    combustion_energy_law.specific_heat_combustion: specific_heat_combustion,
    combustion_energy_law.mass_of_matter: mass_of_fuel
}).rhs

efficiency_equation = efficiency_law.law.subs({
    efficiency_law.active_power: energy_heating_value,
    efficiency_law.full_power: energy_combustion_value,
    efficiency_law.power_factor: efficiency_factor
})
specific_heat_combustion_value = solve(efficiency_equation, specific_heat_combustion, dict=True)[0][specific_heat_combustion]
answer = Eq(specific_heat_combustion, specific_heat_combustion_value)
print(f"Total equation is:\n{print_expression(answer)}")

specific_heat_combustion_j_div_kg = specific_heat_combustion_value.subs({
    efficiency_factor: Quantity(25 * units.percents),
    mass_of_fuel: Quantity(92 * units.grams),
    mass_of_liquid: Quantity(1.8 * units.kilograms),
    temperature_start: to_kelvin_quantity(Celsius(18)),
    temperature_end: to_kelvin_quantity(Celsius(100)),
    specific_heat_heating: Quantity(4200 * units.joules / (units.kilogram * units.kelvin)),
})
answer_value = convert_to(Quantity(specific_heat_combustion_j_div_kg),
                          prefixes.mega * units.joules / units.kilogram)
print(f"Specific heat of combustion is:\n{answer_value} MJ/kg")