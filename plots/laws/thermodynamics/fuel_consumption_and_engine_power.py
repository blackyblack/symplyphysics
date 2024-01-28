#!/usr/bin/env python3

from sympy import symbols, Eq, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.definitions import density_from_mass_volume as density_law
from symplyphysics.laws.kinematic import distance_from_constant_velocity as distance_law
from symplyphysics.laws.electricity import power_factor_from_active_and_full_power as efficiency_law
from symplyphysics.laws.electricity import power_from_energy_time as power_law
from symplyphysics.laws.thermodynamics import energy_from_combustion as combustion_energy_law

# Create efficiency factor values from 15% to 75% in increments of 15%
efficiency_factor_values = [0.15, 0.3, 0.45, 0.6, 0.75]

volume_of_gasoline, velocity_of_car, efficiency_factor, power_of_car, distance = symbols(
    "volume_of_gasoline velocity_of_car efficiency_factor power_of_car distance")
density_of_gasoline, gasoline_specific_heat_combustion = symbols(
    "density_of_gasoline gasoline_specific_heat_combustion")

# Note: fuel_consumption = volume_of_fuel / distance
fuel_consumption = symbols("fuel_consumption")

velocity_equation = distance_law.law.subs({
    distance_law.distance(distance_law.movement_time): distance,
    distance_law.constant_velocity: velocity_of_car,
    distance_law.initial_position: 0
})
time_value = solve(velocity_equation, distance_law.movement_time,
    dict=True)[0][distance_law.movement_time]

power_equation = power_law.law.subs({power_law.time: time_value, power_law.power: power_of_car})
energy_from_power_value = solve(power_equation, power_law.energy, dict=True)[0][power_law.energy]

density_of_gasoline_equation = density_law.definition.subs({
    density_law.density: density_of_gasoline,
    density_law.volume: volume_of_gasoline
})
mass_of_gasoline_value = solve(density_of_gasoline_equation, density_law.mass,
    dict=True)[0][density_law.mass]

amount_heat_value = combustion_energy_law.law.subs({
    combustion_energy_law.specific_heat_combustion: gasoline_specific_heat_combustion,
    combustion_energy_law.mass_of_matter: mass_of_gasoline_value
}).rhs

efficiency_factor_equation = efficiency_law.law.subs({
    efficiency_law.active_power: energy_from_power_value,
    efficiency_law.full_power: amount_heat_value,
    efficiency_law.power_factor: efficiency_factor
})

distance_value = solve(efficiency_factor_equation, distance, dict=True)[0][distance]

fuel_consumption_value = volume_of_gasoline / distance_value
fuel_consumption_equation = Eq(fuel_consumption, fuel_consumption_value)
gasoline_consumption_equation = fuel_consumption_equation.subs({
    velocity_of_car: 15,  # 1 kilometers / hour = (1 / 3.6) meters / second,
    density_of_gasoline: 700,  # kilogram / (meter ** 3),
    gasoline_specific_heat_combustion: 46 * 1E6  # 1 megajoules / kilogram = 1E6 joules / kilogram
})

print(f"Formula is:\n{print_expression(gasoline_consumption_equation)}")

# 1 m^3 / m = 10^(-3) liters / m = (10^(-3) / 10^(-3)) liters / km = 1 liters / km
base_plot = plot(title="Gasoline consumption and engine power",
    xlabel="$Power, watt$",
    ylabel="$V/S, liters/km$",
    backend=MatplotlibBackend,
    legend=True,
    show=False)

# Create plots for every efficient factor in sequence and add plot to base plot
for efficiency_factor_value in efficiency_factor_values:
    gasoline_compustion_value = gasoline_consumption_equation.subs({
        efficiency_factor: efficiency_factor_value
    }).rhs
    subplot = plot(gasoline_compustion_value, (power_of_car, 1_000, 100_000),
        label=r"$\eta_{engine}=" + f"{efficiency_factor_value}$",
        show=False)
    base_plot.append(subplot[0])

base_plot.show()
