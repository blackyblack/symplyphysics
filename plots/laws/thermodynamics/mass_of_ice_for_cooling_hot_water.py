#!/usr/bin/env python3

from sympy import Idx, symbols, Eq, solve, simplify
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression, global_index
from symplyphysics.core.symbols.celsius import to_kelvin, Celsius
from symplyphysics.laws.thermodynamics import thermal_energy_from_heat_capacity_and_temperature as thermal_energy_law
from symplyphysics.laws.thermodynamics import energy_to_melt_from_mass as energy_melting_law
from symplyphysics.definitions import density_from_mass_volume as density_law
from symplyphysics.definitions import specific_heat_capacity as specific_capacity_def
from symplyphysics.laws.thermodynamics import sum_of_heat_transfer_is_zero as thermodinamics_law_1

temperature_of_hot_water_values = [5, 20, 35, 50, 65, 80]

temperature_of_hot_water = symbols("temperature_of_hot_water")
temperature_of_ice = symbols("temperature_of_ice")
temperature_balance = symbols("temperature_balance")

specific_heat_heating_water = symbols("specific_heat_heating_water")
specific_heat_heating_ice = symbols("specific_heat_heating_ice")
specific_heat_melting_ice = symbols("specific_heat_melting_ice")

heat_capacity_expr = solve(
    specific_capacity_def.definition,
    specific_capacity_def.heat_capacity,
)[0]

temperature_melt_ice = symbols("temperature_melt_ice")

mass_of_ice = symbols("mass_of_ice")
mass_of_hot_water = symbols("mass_of_hot_water")

mass_of_all_water = solve(density_law.definition, density_law.symbols.basic.mass,
    dict=True)[0][density_law.symbols.basic.mass]

# the mass of all the water filling the bath consists of the mass of hot water
# that was in the bathroom initially, and the mass of water of melted ice
# mass_all_water = mass_of_hot_water + mass_of_ice
mass_of_all_water_equation = Eq(mass_of_all_water, mass_of_ice + mass_of_hot_water)

energy_cooling_hot_water = thermal_energy_law.law.subs({
    thermal_energy_law.heat_capacity: heat_capacity_expr.subs({
        specific_capacity_def.specific_heat_capacity: specific_heat_heating_water,
        specific_capacity_def.mass: mass_of_hot_water,
    }),
    thermal_energy_law.temperature_origin: temperature_of_hot_water,
    thermal_energy_law.temperature_end: temperature_balance
})

energy_to_heating_ice_equation = thermal_energy_law.law.subs({
    thermal_energy_law.heat_capacity: heat_capacity_expr.subs({
        specific_capacity_def.specific_heat_capacity: specific_heat_heating_ice,
        specific_capacity_def.mass: mass_of_ice,
    }),
    thermal_energy_law.temperature_origin: temperature_of_ice,
    thermal_energy_law.temperature_end: temperature_melt_ice
})

energy_to_melt_ice_equation = energy_melting_law.law.subs({
    energy_melting_law.specific_heat_melting: specific_heat_melting_ice,
    energy_melting_law.symbols.basic.mass: mass_of_ice
})

energy_to_heat_melted_ice_equation = thermal_energy_law.law.subs({
    thermal_energy_law.heat_capacity: heat_capacity_expr.subs({
        specific_capacity_def.specific_heat_capacity: specific_heat_heating_water,
        specific_capacity_def.mass: mass_of_ice,
    }),
    thermal_energy_law.temperature_origin: temperature_melt_ice,
    thermal_energy_law.temperature_end: temperature_balance
})

local_index_ = Idx("local_index_", (1, 4))
thermodinamics_law_1_four_energies = thermodinamics_law_1.law.subs(global_index,
    local_index_).doit()
thermodinamics_law_1_equation = thermodinamics_law_1_four_energies.subs({
    thermodinamics_law_1.amount_energy[1]: energy_cooling_hot_water.rhs,
    thermodinamics_law_1.amount_energy[2]: energy_to_heating_ice_equation.rhs,
    thermodinamics_law_1.amount_energy[3]: energy_to_melt_ice_equation.rhs,
    thermodinamics_law_1.amount_energy[4]: energy_to_heat_melted_ice_equation.rhs,
})

solve_system_equations = solve((thermodinamics_law_1_equation, mass_of_all_water_equation),
    (mass_of_ice, mass_of_hot_water))
mass_of_ice_solve_value = solve_system_equations[mass_of_ice]
mass_of_hot_water_solve_value = solve_system_equations[mass_of_hot_water]

mass_ratio_value = simplify(mass_of_ice_solve_value / mass_of_hot_water_solve_value)
print(f"Total equation:\n{print_expression(Eq(mass_of_ice / mass_of_hot_water, mass_ratio_value))}")

mass_ratio_to_plot = mass_ratio_value.subs({
    temperature_of_ice: to_kelvin(Celsius(-20)),
    specific_heat_heating_water: 4200,  # joules / (kilogram * kelvin)
    specific_heat_heating_ice: 2100,  # joules / (kilogram * kelvin)
    specific_heat_melting_ice: 330_000,  # 1 kilojoules / kilogram = 10^3 joules / kilogram
    temperature_melt_ice: to_kelvin(Celsius(0))
})

base_plot = plot(title="The mass of ice required to cool the hot water to a set temperature",
    xlabel=r"$T_{balance}, K$",
    ylabel=r"$m_{ice} / m_{hot-water}$",
    backend=MatplotlibBackend,
    legend=True,
    show=False)

for temperature_of_hot_water_value in temperature_of_hot_water_values:
    mass_ratio_to_subplot = mass_ratio_to_plot.subs(
        {temperature_of_hot_water: to_kelvin(Celsius(temperature_of_hot_water_value))})

    # Find the upper limit of the temperature scale, so as not to build a graph
    # on the temperature section where there will be a negative mass
    temperature_supremum = solve(mass_ratio_to_subplot.subs({mass_of_ice: 0}),
        temperature_balance,
        dict=True)[0][temperature_balance]

    subplot = plot(mass_ratio_to_subplot,
        (temperature_balance, to_kelvin(Celsius(0)), temperature_supremum),
        label=r"$T_{hot-water}=" + f"{to_kelvin(Celsius(temperature_of_hot_water_value))}" + ", K$",
        show=False)
    base_plot.append(subplot[0])

base_plot.show()
