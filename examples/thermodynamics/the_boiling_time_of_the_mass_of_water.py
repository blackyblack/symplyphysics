#!/usr/bin/env python3

from sympy import solve, Symbol, Eq, Integral
from symplyphysics import print_expression, Quantity, prefixes, units, convert_to
from symplyphysics.core.symbols.celsius import to_kelvin_quantity, Celsius
from symplyphysics.laws.electricity import power_factor_from_active_and_full_power as efficiency_law
from symplyphysics.laws.thermodynamics import energy_from_combustion as combustion_energy_law
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as energy_heating_law
from symplyphysics.laws.thermodynamics import energy_to_vaporization_from_mass as energy_to_vapor_law
from symplyphysics.definitions import mass_flow_rate as mass_rate_law

# Example from https://easyfizika.ru/zadachi/termodinamika/na-zazhzhennuyu-spirtovku-s-kpd-60-postavili-sosud-s-500-g-vody-pri-20-c-cherez-kakoe/
# A vessel with 500 g of water at 20 ° C was placed on a lighted alcohol lamp
# with an efficiency of 60%. After what time will 20 g of water boil off if 4 g of alcohol
# burns in a spirit lamp in a minute?

efficiency = Symbol("efficiency")
mass_of_booling_water = Symbol("mass_of_booling_water")
mass_of_water = Symbol("mass_of_water")
temperature_of_water = Symbol("temperature_of_water")
mass_of_alcohol = Symbol("mass_of_alcohol")
time_of_minute = Symbol("time_of_minute")

mass_rate_alcohol = Symbol("mass_rate_alcohol")

specific_heat_of_combustion_alcohol = Symbol("specific_heat_of_combustion_alcohol")
specific_heat_of_heating_water = Symbol("specific_heat_of_heating_water")
specific_heat_of_vaporization_water = Symbol("specific_heat_of_vaporization_water")
temperature_of_vaporization_water = Symbol("temperature_of_vaporization_water")

time = Symbol("time")

mass_gas_integral = Integral(mass_rate_law.definition,
    (mass_rate_law.time, 0, mass_rate_law.time)).doit()
mass_of_alcohol_in_start_equation = mass_gas_integral.subs({
    mass_rate_law.mass(0): 0,
    mass_rate_law.mass_flow_rate(0): 0,
    mass_rate_law.mass(mass_rate_law.time): mass_of_alcohol,
    mass_rate_law.mass_flow_rate(mass_rate_law.time): mass_rate_alcohol,
    mass_rate_law.time: time_of_minute
}).doit()
mass_flow_rate_in_start_value = solve(mass_of_alcohol_in_start_equation, mass_rate_alcohol,
    dict=True)[0][mass_rate_alcohol]
mass_of_alcohol_value = mass_gas_integral.subs({
    mass_rate_law.mass(0): 0,
    mass_rate_law.mass_flow_rate(0): 0,
    mass_rate_law.mass_flow_rate(mass_rate_law.time): mass_flow_rate_in_start_value,
    mass_rate_law.time: time
}).doit().lhs

energy_from_combustion_alcohol_value = combustion_energy_law.law.subs({
    combustion_energy_law.specific_heat_combustion: specific_heat_of_combustion_alcohol,
    combustion_energy_law.mass_of_matter: mass_of_alcohol_value
}).rhs

energy_for_heating_water_value = energy_heating_law.law.subs({
    energy_heating_law.specific_heat_capacity: specific_heat_of_heating_water,
    energy_heating_law.body_mass: mass_of_water,
    energy_heating_law.temperature_origin: temperature_of_water,
    energy_heating_law.temperature_end: temperature_of_vaporization_water
}).rhs

energy_to_vaporization_water_value = energy_to_vapor_law.law.subs({
    energy_to_vapor_law.mass_of_matter: mass_of_booling_water,
    energy_to_vapor_law.specific_heat_vaporization: specific_heat_of_vaporization_water
}).rhs

# Active work in this example is done to heat water to boiling point and to evaporate. Then
# A_{active} = Q_{heating} + Q_{vaporization}
efficiency_equation = efficiency_law.law.subs({
    efficiency_law.active_power: energy_for_heating_water_value + energy_to_vaporization_water_value,
    efficiency_law.full_power: energy_from_combustion_alcohol_value,
    efficiency_law.power_factor: efficiency
})
time_of_vaporization = solve(efficiency_equation, time,
    dict=True)[0][time]
answer = Eq(time, time_of_vaporization)
print(f"Total equation:\n{print_expression(time_of_vaporization)}")

time_of_vaporization_s = time_of_vaporization.subs({
    efficiency: Quantity(60 * units.percents),
    mass_of_booling_water: Quantity(20 * units.grams),
    mass_of_water: Quantity(500 * units.grams),
    temperature_of_water: to_kelvin_quantity(Celsius(20)),
    mass_of_alcohol: Quantity(4 * units.grams),
    time_of_minute: Quantity(1 * units.minutes),

    specific_heat_of_combustion_alcohol: Quantity(29 * prefixes.mega * units.joules / units.kilogram),
    specific_heat_of_heating_water: Quantity(4200 * units.joules / (units.kilogram * units.kelvin)),
    specific_heat_of_vaporization_water: Quantity(2.26 * prefixes.mega * units.joules / units.kilogram),
    temperature_of_vaporization_water: to_kelvin_quantity(Celsius(100)),
})
print(f"Time for vaporization water is: {convert_to(Quantity(time_of_vaporization_s), units.seconds)} s")
