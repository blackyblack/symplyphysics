#!/usr/bin/env python3
"""
Plot the fuel consumption by a car engine as a function of the engine power for different values
of the engine's efficiency factor.

**Fuel consumption** is defined as the volume of the fuel consumed per unit distance traveled.
"""

from sympy import symbols, Eq, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, units, Quantity
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.classical_mechanics.fundamentals import density_from_mass_volume as density_law
from symplyphysics.classical_mechanics.kinematics.translational_motion import position_via_constant_speed_and_time as distance_law
from symplyphysics.electromagnetism.circuits.alternating_current import power_factor_is_real_power_over_apparent_power as efficiency_law
from symplyphysics.classical_mechanics.fundamentals import energy_via_constant_power_and_time as energy_law
from symplyphysics.thermodynamics.phase_transitions.latent_heat import heat_of_combustion_via_mass as combustion_energy_law

efficiency_factor_values = [0.15, 0.3, 0.45, 0.6, 0.75]

volume_of_gasoline, velocity_of_car, efficiency_factor, power_of_car, distance = symbols(
    "volume_of_gasoline velocity_of_car efficiency_factor power_of_car distance")
density_of_gasoline, gasoline_specific_heat_combustion = symbols(
    "density_of_gasoline gasoline_specific_heat_combustion")

fuel_consumption = symbols("fuel_consumption")

velocity_equation = distance_law.law.subs({
    distance_law.final_position: distance,
    distance_law.speed: velocity_of_car,
    distance_law.initial_position: 0
})
time_value = solve(velocity_equation, distance_law.time)[0]

energy_from_power_value = energy_law.law.rhs.subs({
    energy_law.time: time_value,
    energy_law.power: power_of_car,
})

density_of_gasoline_equation = density_law.law.subs({
    density_law.density: density_of_gasoline,
    density_law.volume: volume_of_gasoline
})
mass_of_gasoline_value = solve(density_of_gasoline_equation, density_law.mass)[0]

amount_heat_value = combustion_energy_law.law.subs({
    combustion_energy_law.specific_heat_of_combustion: gasoline_specific_heat_combustion,
    combustion_energy_law.mass: mass_of_gasoline_value
}).rhs

efficiency_factor_equation = efficiency_law.law.subs({
    efficiency_law.real_power: energy_from_power_value,
    efficiency_law.apparent_power: amount_heat_value,
    efficiency_law.power_factor: efficiency_factor
})

distance_value = solve(efficiency_factor_equation, distance)[0]

fuel_consumption_expr = volume_of_gasoline / distance_value
fuel_consumption_equation = Eq(fuel_consumption, fuel_consumption_expr)
gasoline_consumption_equation = fuel_consumption_equation.subs({
    velocity_of_car: 54 * units.kilometer / units.hour,
    density_of_gasoline: 700 * units.kilogram / units.meter**3,
    gasoline_specific_heat_combustion: Quantity(46 * units.mega * units.joule / units.kilogram),
})

consumption_in_l_per_km = symbols("k")
consumption_eqn = Eq(fuel_consumption, consumption_in_l_per_km * units.liter / units.kilometer)

power_in_kw = symbols("p")
power_eqn = Eq(power_of_car, power_in_kw * Quantity(units.kilo * units.watt))

consumption_expr = solve(
    (gasoline_consumption_equation, consumption_eqn, power_eqn),
    (fuel_consumption, consumption_in_l_per_km, power_of_car),
    dict=True,
)[0][consumption_in_l_per_km]
consumption_expr = evaluate_expression(consumption_expr)

print(
    "Fuel consumption (in L/km) as a function of power `p` (in W):",
    print_expression(consumption_expr),
    sep="\n",
)

base_plot = plot(
    title="Gasoline consumption and engine power",
    xlabel="power, kW",
    ylabel="fuel consumption, L/km",
    backend=MatplotlibBackend,
    legend=True,
    show=False,
)

# Create plots for every efficient factor in sequence and add plot to base plot
for efficiency_factor_value in efficiency_factor_values:
    consumption_expr_subs = consumption_expr.subs(efficiency_factor, efficiency_factor_value)
    subplot = plot(
        consumption_expr_subs,
        (power_in_kw, 1, 100),
        label=rf"$\eta_\text{{engine}} = {efficiency_factor_value}$",
        show=False,
    )
    base_plot.extend(subplot)

base_plot.show()
