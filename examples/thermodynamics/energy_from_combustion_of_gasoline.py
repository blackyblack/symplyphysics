#!/usr/bin/env python3

from sympy import solve, Symbol, Eq
from symplyphysics import print_expression, Quantity, prefixes, units, convert_to
from symplyphysics.definitions import density_from_mass_volume as density_law
from symplyphysics.laws.kinematics import position_via_constant_speed_and_time as distance_law
from symplyphysics.laws.electricity import power_factor_from_active_and_full_power as efficiency_law
from symplyphysics.laws.electricity import energy_via_constant_power_and_time as energy_law
from symplyphysics.laws.thermodynamics import heat_of_combustion_via_mass as combustion_energy_law

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

velocity_equation = distance_law.law.subs({
    distance_law.final_position: distance,
    distance_law.speed: velocity_of_car,
    distance_law.initial_position: 0
})
time_value = solve(velocity_equation, distance_law.time,
    dict=True)[0][distance_law.time]

# power_equation = power_law.law.subs({power_law.time: time_value, power_law.power: power_of_car})
energy_from_power_value = energy_law.law.rhs.subs({
    energy_law.time: time_value,
    energy_law.power: power_of_car,
})

density_of_gasoline_equation = density_law.definition.subs({
    density_law.density: density_of_gasoline,
    density_law.volume: volume_of_gasoline
})
mass_of_gasoline_value = solve(density_of_gasoline_equation, density_law.mass,
    dict=True)[0][density_law.mass]

amount_heat_value = combustion_energy_law.law.subs({
    combustion_energy_law.specific_heat_of_combustion: gasoline_specific_heat_combustion,
    combustion_energy_law.mass: mass_of_gasoline_value
}).rhs

efficiency_factor_equation = efficiency_law.law.subs({
    efficiency_law.active_power: energy_from_power_value,
    efficiency_law.full_power: amount_heat_value,
    efficiency_law.power_factor: efficiency_factor
})
distance_value = solve(efficiency_factor_equation, distance, dict=True)[0][distance]
answer = Eq(distance, distance_value)
print(f"Total equation:\n{print_expression(answer)}")

distance_m = distance_value.subs({
    volume_of_gasoline: Quantity(40 * units.liters),
    velocity_of_car: Quantity(54 * units.kilometers / units.hour),
    power_of_car: Quantity(46 * prefixes.kilo * units.watts),
    efficiency_factor: Quantity(25 * units.percents),
    density_of_gasoline: Quantity(700 * units.kilogram / (units.meter**3)),
    gasoline_specific_heat_combustion: Quantity(46 * prefixes.mega * units.joules / units.kilogram)
})
print(f"Distance is: {convert_to(Quantity(distance_m), units.kilometers).evalf(4)} km")
