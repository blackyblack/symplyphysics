#!/usr/bin/env python3

from sympy import symbols, Eq, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as energy_heating_law
from symplyphysics.laws.thermodynamics import energy_to_melt_from_mass as energy_melting_law
from symplyphysics.definitions import density_from_mass_volume as density_law
from symplyphysics.laws.thermodynamics import sum_of_heat_transfer_is_zero as thermodinamics_law_1

bath_volumes = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]

volume_of_bath = symbols("volume_of_bath")
temperature_of_hot_water = symbols("temperature_of_hot_water")
temperature_of_ice = symbols("temperature_of_ice")
temperature_balance = symbols("temperature_balance")

density_of_water = symbols("density_of_water")
specific_heat_heating_water = symbols("specific_heat_heating_water")
specific_heat_heating_ice = symbols("specific_heat_heating_ice")
specific_heat_melting_ice = symbols("specific_heat_melting_ice")
temperature_melt_ice = symbols("temperature_melt_ice")

mass_of_ice = symbols("mass_of_ice")
mass_of_hot_water = symbols("mass_of_hot_water")

density_of_water_equation = density_law.definition.subs({
    density_law.volume: volume_of_bath,
    density_law.density: density_of_water
})
mass_of_all_water = solve(density_of_water_equation, density_law.mass,
    dict=True)[0][density_law.mass]

# the mass of all the water filling the bath consists of the mass of hot water
# that was in the bathroom initially, and the mass of water of melted ice
# mass_all_water = mass_of_hot_water + mass_of_ice
mass_of_all_water_equation = Eq(mass_of_all_water, mass_of_ice + mass_of_hot_water)
mass_of_hot_water_value = solve(mass_of_all_water_equation, mass_of_hot_water,
    dict=True)[0][mass_of_hot_water]

energy_cooling_hot_water = energy_heating_law.law.subs({
    energy_heating_law.specific_heat_capacity: specific_heat_heating_water,
    energy_heating_law.body_mass: mass_of_hot_water_value,
    energy_heating_law.temperature_origin: temperature_of_hot_water,
    energy_heating_law.temperature_end: temperature_balance
})

energy_to_heating_ice_equation = energy_heating_law.law.subs({
    energy_heating_law.specific_heat_capacity: specific_heat_heating_ice,
    energy_heating_law.body_mass: mass_of_ice,
    energy_heating_law.temperature_origin: temperature_of_ice,
    energy_heating_law.temperature_end: temperature_melt_ice
})

energy_to_melt_ice_equation = energy_melting_law.law.subs({
    energy_melting_law.specific_heat_melting: specific_heat_melting_ice,
    energy_melting_law.mass_of_matter: mass_of_ice
})

energy_to_heat_melted_ice_equation = energy_heating_law.law.subs({
    energy_heating_law.specific_heat_capacity: specific_heat_heating_water,
    energy_heating_law.body_mass: mass_of_ice,
    energy_heating_law.temperature_origin: temperature_melt_ice,
    energy_heating_law.temperature_end: temperature_balance
})

thermodinamics_law_1_equation = thermodinamics_law_1.law.subs({
    thermodinamics_law_1.amounts_energy:
    (energy_cooling_hot_water.rhs, energy_to_heating_ice_equation.rhs,
    energy_to_melt_ice_equation.rhs, energy_to_heat_melted_ice_equation.rhs)
}).doit()

mass_of_ice_value = solve(thermodinamics_law_1_equation, mass_of_ice, dict=True)[0][mass_of_ice]
mass_of_ice_equation = Eq(mass_of_ice, mass_of_ice_value)
print(f"Formula is:\n{print_expression(mass_of_ice_equation)}")

mass_of_ice_equation_to_plot = mass_of_ice_equation.subs({
    temperature_of_hot_water: 273.15 + 80,      # 80 degree Celsius = 273.15 + 80 kelvins
    temperature_of_ice: 273.15 - 20,        # -20 degree Celsius = 273.15 - 20 kelvins
    density_of_water: 1000,     # kilograms / (meter^3)
    specific_heat_heating_water: 4200,      # joules / (kilogram * kelvin)
    specific_heat_heating_ice: 2100,        # joules / (kilogram * kelvin)
    specific_heat_melting_ice: 330_000,     # 1 kilojoules / kilogram = 10^3 joules / kilogram
    temperature_melt_ice: 273.15 + 0
})

base_plot = plot(title="The mass of ice required to cool the hot water to a set temperature",
    xlabel=r"$T_\text{balance}, K$",
    ylabel=r"$m_\text{ice}, kg$",
    backend=MatplotlibBackend,
    legend=True,
    show=False)

for bath_volume in bath_volumes:
    mass_of_ice_equation_to_subplot = mass_of_ice_equation_to_plot.subs({
        volume_of_bath: bath_volume
    })

    # Find the upper limit of the temperature scale, so as not to build a graph
    # on the temperature section where there will be a negative mass
    temperature_supremum = solve(mass_of_ice_equation_to_subplot.subs({mass_of_ice: 0}), temperature_balance,
        dict=True)[0][temperature_balance]

    subplot = plot(mass_of_ice_equation_to_subplot.rhs, (temperature_balance, 280, temperature_supremum),
        label=r"$V_\text{bath}=" + f"{bath_volume}" + "\ m^3$",
        show=False)
    base_plot.append(subplot[0])

base_plot.show()
