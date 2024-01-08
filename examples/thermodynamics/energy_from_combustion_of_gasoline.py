from sympy import solve, Symbol, Eq
from symplyphysics import print_expression
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

velocity_equation = distance_law.law.subs({
    distance_law.distance(distance_law.movement_time): distance,
    distance_law.constant_velocity: velocity_of_car,
    distance_law.initial_position: 0
})
time_equation = solve(velocity_equation, distance_law.movement_time, dict=True)[0][distance_law.movement_time]

power_equation = power_law.law.subs({
    power_law.time: time_equation,
    power_law.power: power_of_car
})
energy_from_power_equation = solve(power_equation, power_law.energy, dict=True)[0][power_law.energy]

density_of_gasoline_equation = density_law.definition.subs({
    density_law.density: density_of_gasoline,
    density_law.volume: volume_of_gasoline
})
mass_of_gasoline_equation = solve(density_of_gasoline_equation, density_law.mass, dict=True)[0][density_law.mass]

amount_heat_equation = combustion_energy_law.law.subs({
    combustion_energy_law.specific_heat_combustion: gasoline_specific_heat_combustion,
    combustion_energy_law.mass_of_matter: mass_of_gasoline_equation
}).rhs

efficiency_factor_equation = efficiency_law.law.subs({
    efficiency_law.active_power: energy_from_power_equation,
    efficiency_law.full_power: amount_heat_equation,
    efficiency_law.power_factor: efficiency_factor
})
distance_equation = solve(efficiency_factor_equation, distance, dict=True)[0][distance]
answer = Eq(distance, distance_equation)
print(f"Total equation:\n{print_expression(answer)}")

distance_m = distance_equation.subs({
    volume_of_gasoline: 40 * 1E-3,
    velocity_of_car: 15,
    power_of_car: 46 * 1E3,
    efficiency_factor: 0.25,
    density_of_gasoline: 700,
    gasoline_specific_heat_combustion: 46 * 1E6

})
print(f"Distance is: {distance_m} m")
