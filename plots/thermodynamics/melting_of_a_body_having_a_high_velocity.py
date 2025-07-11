#!/usr/bin/env python3

from sympy import symbols, Eq, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.core.symbols.celsius import to_kelvin, Celsius
from symplyphysics.laws.thermodynamics import heat_is_heat_capacity_times_temperature_change as thermal_energy_law
from symplyphysics.laws.quantities import quantity_is_specific_quantity_times_mass as specific_qty_law
from symplyphysics.laws.thermodynamics import latent_heat_of_fusion_via_mass as energy_melting_law
from symplyphysics.laws.conservation import initial_mechanical_energy_equals_final_mechanical_energy as energy_conservation_law
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as kinetic_energy_law

matter_parameters: dict[str, dict[str, float]] = {
    "Fe": {
        "specific_heat_heating": 460,
        "specific_heat_melting": 270_000,
        "temperature_of_melting": 1_400,
    },
    "W": {
        "specific_heat_heating": 134,
        "specific_heat_melting": 184_000,
        "temperature_of_melting": 3_422,
    },
    "Al": {
        "specific_heat_heating": 897,
        "specific_heat_melting": 393_000,
        "temperature_of_melting": 660,
    },
    "Si": {     # silicon/silicium
        "specific_heat_heating": 714,
        "specific_heat_melting": 1_409_000,
        "temperature_of_melting": 1_415,
    },
    "Cu": {
        "specific_heat_heating": 383,
        "specific_heat_melting": 213_000,
        "temperature_of_melting": 1_083.4,
    },
    "Au": {
        "specific_heat_heating": 129,
        "specific_heat_melting": 67_000,
        "temperature_of_melting": 1_063,
    },
    "Zn": {
        "specific_heat_heating": 400,
        "specific_heat_melting": 117_000,
        "temperature_of_melting": 415,
    }
}

velocity_of_meteorite = symbols("velocity_of_meteorite")
temperature_of_meteorite = symbols("temperature_of_meteorite")
koefficient_of_transfer_to_kinetic = symbols("koefficient_of_transfer_to_kinetic")

specific_heat_heating_meteorite = symbols("specific_heat_heating_meteorite")
specific_heat_melting_meteorite = symbols("specific_heat_melting_meteorite")
temperature_of_meteorite_melting = symbols("temperature_of_meteorite_melting")

mass_of_meteorite = symbols("mass_of_meteorite")
mass_of_melting_meteorite = symbols("mass_of_melting_meteorite")
koefficient_of_melting_meteorite = symbols("koefficient_of_melting_meteorite")

meteorite_heat_capacity = specific_qty_law.law.rhs.subs({
    specific_qty_law.specific_quantity: specific_heat_heating_meteorite,
    specific_qty_law.mass: mass_of_meteorite,
})

energy_to_heating_meteorite = thermal_energy_law.law.subs({
    thermal_energy_law.heat_capacity:
        meteorite_heat_capacity,
    thermal_energy_law.temperature_change:
    temperature_of_meteorite_melting - temperature_of_meteorite,
}).rhs

energy_to_meteorite_melting = energy_melting_law.law.subs({
    energy_melting_law.mass: mass_of_melting_meteorite,
    energy_melting_law.specific_heat_of_fusion: specific_heat_melting_meteorite
}).rhs

kinetic_energy_of_meteorite = kinetic_energy_law.law.subs({
    kinetic_energy_law.mass: mass_of_meteorite,
    kinetic_energy_law.speed: velocity_of_meteorite
}).rhs

energy_conservation_equation = energy_conservation_law.law.subs({
    energy_conservation_law.mechanical_energy(energy_conservation_law.initial_time):
    koefficient_of_transfer_to_kinetic * kinetic_energy_of_meteorite,
    energy_conservation_law.mechanical_energy(energy_conservation_law.final_time):
    energy_to_heating_meteorite + energy_to_meteorite_melting
})
mass_of_melting_value = solve(energy_conservation_equation, mass_of_melting_meteorite,
    dict=True)[0][mass_of_melting_meteorite]
koefficient_of_melting_meteorite_value = mass_of_melting_value / mass_of_meteorite
answer = Eq(koefficient_of_melting_meteorite, koefficient_of_melting_meteorite_value)
print(f"Total equation:\n{print_expression(answer)}")

koefficient_of_melting_meteorite_to_plots = answer.subs({
    koefficient_of_transfer_to_kinetic: 0.80,  # percents
    temperature_of_meteorite: 300  # kelvins
})

base_plot = plot(title="The proportion of a molten meteorite depending on its velocity",
    xlabel=r"$v, \text{m}/\text{s}$",
    ylabel=r"$\text{molten ratio}$",
    backend=MatplotlibBackend,
    legend=True,
    show=False)

MASSES_RATIO_MAXIMUM = 1
MASSES_RATIO_MINIMUM = 0
for matter, parameters in matter_parameters.items():
    masses_ratio_to_subplot = koefficient_of_melting_meteorite_to_plots.subs({
        specific_heat_heating_meteorite:
        parameters["specific_heat_heating"],  # joules / (kilogram * kelvin)
        specific_heat_melting_meteorite: parameters["specific_heat_melting"],  # joules / kilogram
        temperature_of_meteorite_melting: to_kelvin(Celsius(parameters["temperature_of_melting"])),
    })
    # First solve is negative value. Ignore it.
    masses_ratio_maximum_equation = masses_ratio_to_subplot.subs(
        {koefficient_of_melting_meteorite: MASSES_RATIO_MAXIMUM})
    velocity_maximum = solve(masses_ratio_maximum_equation, velocity_of_meteorite,
        dict=True)[-1][velocity_of_meteorite]
    masses_ratio_minimum_equation = masses_ratio_to_subplot.subs(
        {koefficient_of_melting_meteorite: MASSES_RATIO_MINIMUM})
    velocity_minimum = solve(masses_ratio_minimum_equation, velocity_of_meteorite,
        dict=True)[-1][velocity_of_meteorite]

    subplot = plot(masses_ratio_to_subplot.rhs,
        (velocity_of_meteorite, velocity_minimum, velocity_maximum),
        label=f"Material of meteorite is {matter}",
        show=False)
    base_plot.append(subplot[0])

base_plot.show()
