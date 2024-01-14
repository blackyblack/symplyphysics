#!/usr/bin/env python3

from sympy import dsolve, solve, Symbol, Eq, pi
from symplyphysics import print_expression, Quantity, prefixes, units, convert_to
from symplyphysics.core.symbols.celsius import to_kelvin_quantity, Celsius
from symplyphysics.laws.electricity import power_factor_from_active_and_full_power as efficiency_law
from symplyphysics.laws.thermodynamics import energy_from_combustion as combustion_energy_law
from symplyphysics.definitions import density_from_mass_volume as density_law
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as energy_heating_law
from symplyphysics.laws.kinematic import distance_from_constant_velocity as velocity_law
from symplyphysics.laws.thermodynamics import pressure_from_temperature_and_volume as klayperon_law
from symplyphysics.laws.chemistry import atomic_weight_from_mass_mole_count as mole_count_law
from symplyphysics.definitions import mass_flow_rate as mass_rate_law

# Example from https://easyfizika.ru/zadachi/termodinamika/gazovaya-nagrevatelnaya-kolonka-potreblyaet-1-8-m3-metana-ch4-v-chas-najti-temperaturu/
# The gas heating column consumes 1.8 m^3 of methane (CH4) per hour. Find
# the temperature of the water flowing out, if its jet has a velocity of 0.5 m / s
# and a diameter of 1 cm, the initial temperature of water and gas is 11 ° C,
# the calorific value of methane is 55 kJ / g. The gas in the pipe is at a pressure of 120 kPa.
# The efficiency of the heater is 60%.

volume_of_gas = Symbol("volume_of_gas")
diameter_of_pipe = Symbol("diameter_of_pipe")
time_of_hour = Symbol("time_of_hour")
velocity_of_water = Symbol("velocity_of_water")
pressure_in_gas_heater = Symbol("pressure_in_gas_heater")
efficiency_of_gas_heater = Symbol("efficiency_of_gas_heater")
temperature_start = Symbol("temperature_start")

mass_flow_rate = Symbol("mass_flow_rate")

molar_mass_of_metan = Symbol("molar_mass_of_metan")
specific_heat_of_combustion_metan = Symbol("specific_heat_of_combustion_metan")
specific_heat_of_heating_water = Symbol("specific_heat_of_heating_water")
density_of_water = Symbol("density_of_water")

temperature_water = Symbol("temperature_water")

distance_value = velocity_law.law.subs({
    velocity_law.constant_velocity: velocity_of_water,
    velocity_law.initial_position: 0
}).rhs

# We assume that the cross-section of the pipe through which the water flows has the shape of a circle
# S = pi * d^2 / 4
pipe_cross_sectional_area = pi * diameter_of_pipe ** 2 / 4
# And we assume that the pipe goes straight to the water heating site.
# Accordingly, the volume of the pipe that the water passes over a certain period of time has the shape of a cylinder
# V = S * l(t)
volume_value = pipe_cross_sectional_area * distance_value

density_of_water_equation = density_law.definition.subs({
    density_law.density: density_of_water,
    density_law.volume: volume_value
})
mass_of_water_value = solve(density_of_water_equation, density_law.mass,
    dict=True)[0][density_law.mass]

energy_to_heating_water_value = energy_heating_law.law.subs({
    energy_heating_law.body_mass: mass_of_water_value,
    energy_heating_law.specific_heat_capacity: specific_heat_of_heating_water,
    energy_heating_law.temperature_origin: temperature_start,
    energy_heating_law.temperature_end: temperature_water
}).rhs

mass_of_gas_equation = mole_count_law.law.subs({
    mole_count_law.atomic_weight: molar_mass_of_metan
})
mole_count_value = solve(mass_of_gas_equation, mole_count_law.mole_count,
    dict=True)[0][mole_count_law.mole_count]

state_equation = klayperon_law.law.subs({
    klayperon_law.volume: volume_of_gas,
    klayperon_law.pressure: pressure_in_gas_heater,
    klayperon_law.temperature: temperature_start,
    klayperon_law.mole_count: mole_count_value
})
mass_of_gas_in_state_value = solve(state_equation, mole_count_law.substance_mass,
    dict=True)[0][mole_count_law.substance_mass]

mass_gas_dsolved = dsolve(mass_rate_law.definition, mass_rate_law.mass(mass_rate_law.time))
# C1 is initial mass of consumed gas
# HACK: we should tell dsolve() that mass flow rate is constant, but solve() does not work properly with constant
#       value. So we substitute it with constant and revert to normal mass flow rate after dsolve()
mass_flow_constant_rate = Symbol("mass_flow_rate_const", constant=True)
mass_gas_eq = mass_gas_dsolved.subs({"C1": 0, mass_rate_law.mass_flow_rate(mass_rate_law.time): mass_flow_constant_rate}).doit()
mass_gas_eq = mass_gas_eq.subs(mass_flow_constant_rate, mass_flow_rate)

mass_of_gas_in_start_equation = mass_gas_eq.subs({
    mass_rate_law.mass(mass_rate_law.time): mass_of_gas_in_state_value,
    mass_rate_law.time: time_of_hour
})
mass_flow_rate_in_start_value = solve(mass_of_gas_in_start_equation, mass_flow_rate,
    dict=True)[0][mass_flow_rate]
mass_of_gas_value = solve(mass_gas_eq.subs({
    mass_flow_rate: mass_flow_rate_in_start_value,
    mass_rate_law.time: velocity_law.movement_time
}), mass_rate_law.mass(velocity_law.movement_time))[0]

energy_from_combustion_of_metan_value = combustion_energy_law.law.subs({
    combustion_energy_law.specific_heat_combustion: specific_heat_of_combustion_metan,
    combustion_energy_law.mass_of_matter: mass_of_gas_value
}).rhs

efficiency_equation = efficiency_law.law.subs({
    efficiency_law.active_power: energy_to_heating_water_value,
    efficiency_law.full_power: energy_from_combustion_of_metan_value,
    efficiency_law.power_factor: efficiency_of_gas_heater
})

temperature_water_value = solve(efficiency_equation, temperature_water,
    dict=True)[0][temperature_water]
answer = Eq(temperature_water, temperature_water_value)
print(f"Total equation:\n{print_expression(answer)}")

temperature_water_k = temperature_water_value.subs({
    volume_of_gas: Quantity(1.8 * units.meters ** 3),
    diameter_of_pipe: Quantity(1 * prefixes.centi * units.meters),
    time_of_hour: Quantity(1 * units.hours),
    velocity_of_water: Quantity(0.5 * units.meters / units.second),
    pressure_in_gas_heater: Quantity(120 * prefixes.kilo * units.pascals),
    efficiency_of_gas_heater: Quantity(60 * units.percents),
    temperature_start: to_kelvin_quantity(Celsius(11)),
    molar_mass_of_metan: Quantity(16 * units.grams / units.mole),
    specific_heat_of_combustion_metan: Quantity(55 * prefixes.mega * units.joules / units.kilogram),
    specific_heat_of_heating_water: Quantity(4200 * units.joules / (units.kilogram * units.kelvin)),
    density_of_water: Quantity(1_000 * units.kilograms / (units.meter ** 3))
})
print(f"Temperature of water in out is: {convert_to(Quantity(temperature_water_k), units.kelvins).evalf(5)} K")
