#!/usr/bin/env python3

from sympy import symbols, Eq, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.core.symbols.celsius import to_kelvin, Celsius
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as energy_heating_law
from symplyphysics.laws.thermodynamics import energy_to_melt_from_mass as energy_melting_law
from symplyphysics.laws.conservation import mechanical_energy_after_equals_to_mechanical_energy_before as energy_conservation_law
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy_law

matter_parameters = {
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

energy_to_heating_meteorite = energy_heating_law.law.subs({
    energy_heating_law.specific_heat_capacity: specific_heat_heating_meteorite,
    energy_heating_law.body_mass: mass_of_meteorite,
    energy_heating_law.temperature_origin: temperature_of_meteorite,
    energy_heating_law.temperature_end: temperature_of_meteorite_melting
}).rhs

energy_to_meteorite_melting = energy_melting_law.law.subs({
    energy_melting_law.mass_of_matter: mass_of_melting_meteorite,
    energy_melting_law.specific_heat_melting: specific_heat_melting_meteorite
}).rhs

kinetic_energy_of_meteorite = kinetic_energy_law.law.subs({
    kinetic_energy_law.body_mass: mass_of_meteorite,
    kinetic_energy_law.body_velocity: velocity_of_meteorite
}).rhs

energy_conservation_equation = energy_conservation_law.law.subs({
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_before): koefficient_of_transfer_to_kinetic * kinetic_energy_of_meteorite,
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_after): energy_to_heating_meteorite + energy_to_meteorite_melting
})
mass_of_melting_value = solve(energy_conservation_equation, mass_of_melting_meteorite,
    dict=True)[0][mass_of_melting_meteorite]
koefficient_of_melting_meteorite_value = mass_of_melting_value / mass_of_meteorite
answer = Eq(koefficient_of_melting_meteorite, koefficient_of_melting_meteorite_value)
print(f"Total equation:\n{print_expression(answer)}")

koefficient_of_melting_meteorite_to_plots = answer.subs({
    koefficient_of_transfer_to_kinetic: 0.80,   # percents
    temperature_of_meteorite: 300   # kelvins
})

base_plot = plot(title="The proportion of a molten meteorite depending on its velocity",
    xlabel=r"$v, m/s$",
    ylabel=r"$molten ratio$",
    backend=MatplotlibBackend,
    legend=True,
    show=False)

MASSES_RATIO_MAXIMUM = 1
MASSES_RATIO_MINIMUM = 0
for matter, parameters in matter_parametrs.items():
    masses_ratio_to_subplot = koefficient_of_melting_meteorite_to_plots.subs({
        specific_heat_heating_meteorite: dict(parameters)["specific_heat_heating"],   # joules / (kilogram * kelvin)
        specific_heat_melting_meteorite: dict(parameters)["specific_heat_melting"],   # joules / kilogram
        temperature_of_meteorite_melting: to_kelvin(Celsius(dict(parameters)["temperature_of_melting"]))
    })
    # First solve is negative value. Ignore it.
    masses_ratio_supremum_equation = masses_ratio_to_subplot.subs({
        koefficient_of_melting_meteorite: MASSES_RATIO_MAXIMUM
    })
    velocity_maximum = solve(masses_ratio_supremum_equation, velocity_of_meteorite,
        dict=True)[-1][velocity_of_meteorite]
    masses_ratio_minimum_equation = masses_ratio_to_subplot.subs({
        koefficient_of_melting_meteorite: MASSES_RATIO_MINIMUM
    })
    velocity_minimum = solve(masses_ratio_minimum_equation, velocity_of_meteorite,
        dict=True)[-1][velocity_of_meteorite]

    subplot = plot(masses_ratio_to_subplot.rhs, (velocity_of_meteorite, velocity_minimum, velocity_maximum),
                   label=f"Material of meteorite is {matter}",
                   show=False)
    base_plot.append(subplot[0])

base_plot.show()
