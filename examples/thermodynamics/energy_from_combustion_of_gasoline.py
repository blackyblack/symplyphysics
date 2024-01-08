#!/usr/bin/env python3

from sympy import solve, Symbol, Eq
from symplyphysics import print_expression, Quantity, prefixes, units, \
    convert_to, dimensionless
from symplyphysics.definitions import density_from_mass_volume as density_law
from symplyphysics.laws.kinematic import distance_from_constant_velocity as distance_law
from symplyphysics.laws.electricity import power_factor_from_active_and_full_power as efficiency_law
from symplyphysics.laws.electricity import power_from_energy_time as power_law
from symplyphysics.laws.thermodynamics import energy_from_combustion as combustion_energy_law

# Example from https://easyfizika.ru/zadachi/termodinamika/na-skolko-kilometrov-puti-hvatit-40-l-benzina-avtomobilyu-dvizhushhemusya-so-skorostyu/
# How many kilometers will 40 liters of gasoline be enough for a car moving at a speed of 54 km/h?
# The power developed by the car engine is 46 kW, and the efficiency is 25%.

volume_of_gasoline = Symbol("volume_of_gasoline")
velocity_of_car = Symbol("velocity_of_car")
efficiency_factor = Symbol("efficiency_factor")
power_of_car = Symbol("power_of_car")

density_of_gasoline = Symbol("density_of_gasoline")
gasoline_specific_heat_combustion = Symbol("gasoline_specific_heat_combustion")

distance = Symbol("distance")

velocity_value = distance_law.law.subs({
    distance_law.distance(distance_law.movement_time): distance,
    distance_law.constant_velocity: velocity_of_car,
    distance_law.initial_position: 0
})
time_value = solve(velocity_value, distance_law.movement_time, dict=True)[0][distance_law.movement_time]

power_value = power_law.law.subs({
    power_law.time: time_value,
    power_law.power: power_of_car
})
energy_from_power_value = solve(power_value, power_law.energy, dict=True)[0][power_law.energy]

density_of_gasoline_value = density_law.definition.subs({
    density_law.density: density_of_gasoline,
    density_law.volume: volume_of_gasoline
})
mass_of_gasoline_value = solve(density_of_gasoline_value, density_law.mass, dict=True)[0][density_law.mass]

amount_heat_value = combustion_energy_law.law.subs({
    combustion_energy_law.specific_heat_combustion: gasoline_specific_heat_combustion,
    combustion_energy_law.mass_of_matter: mass_of_gasoline_value
}).rhs

efficiency_factor_value = efficiency_law.law.subs({
    efficiency_law.active_power: energy_from_power_value,
    efficiency_law.full_power: amount_heat_value,
    efficiency_law.power_factor: efficiency_factor
})
distance_value = solve(efficiency_factor_value, distance, dict=True)[0][distance]
answer = Eq(distance, distance_value)
print(f"Total equation:\n{print_expression(answer)}")

distance_m = distance_value.subs({
    volume_of_gasoline: convert_to(Quantity(40 * units.liters), units.meters ** 3),
    velocity_of_car: convert_to(Quantity(54 * units.kilometers / units.hour), units.meters / units.second),
    power_of_car: convert_to(Quantity(46 * prefixes.kilo * units.watts), units.watts),
    efficiency_factor: convert_to(Quantity(25 * units.percents), dimensionless),
    density_of_gasoline: convert_to(Quantity(700 * units.kilogram / (units.meter ** 3)), units.kilogram / (units.meter ** 3)),
    gasoline_specific_heat_combustion: convert_to(Quantity(46 * prefixes.mega * units.joules / units.kilogram), units.joules / units.kilogram)
})
print(f"Distance is: {distance_m} m")
