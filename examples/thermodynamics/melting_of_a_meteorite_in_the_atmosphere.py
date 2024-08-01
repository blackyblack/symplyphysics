#!/usr/bin/env python3

from sympy import solve, Symbol, Eq
from symplyphysics import print_expression, Quantity, prefixes, units, convert_to
from symplyphysics.core.symbols.celsius import to_kelvin_quantity, Celsius
from symplyphysics.laws.thermodynamics import (
    latent_heat_of_fusion_via_mass as energy_melting_law,
    thermal_energy_from_heat_capacity_and_temperature as thermal_energy_law,
)
from symplyphysics.laws.quantities import quantity_is_specific_quantity_times_mass as specific_qty_law
from symplyphysics.laws.conservation import (
    mechanical_energy_after_equals_to_mechanical_energy_before as energy_conservation_law,)
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy_law

# Example from https://easyfizika.ru/zadachi/termodinamika/zheleznyj-meteorit-vletaet-v-atmosferu-zemli-so-skorostyu-1-5x10-3-m-s/
# An iron meteorite flies into the Earth's atmosphere at a speed of 1.5 km/s,
# having a temperature of 300 K. If 80% of the kinetic energy of a meteorite passes
# into its internal energy when moving in the atmosphere,
# then what part of the meteorite's mass will melt.

velocity_of_meteorite = Symbol("velocity_of_meteorite")
temperature_of_meteorite = Symbol("temperature_of_meteorite")
koefficient_of_transfer_to_kinetic = Symbol("koefficient_of_transfer_to_kinetic")

specific_heat_heating_meteorite = Symbol("specific_heat_heating_meteorite")
specific_heat_melting_meteorite = Symbol("specific_heat_melting_meteorite")
temperature_of_meteorite_melting = Symbol("temperature_of_meteorite_melting")

mass_of_meteorite = Symbol("mass_of_meteorite")
mass_of_melting_meteorite = Symbol("mass_of_melting_meteorite")
koefficient_of_melting_meteorite = Symbol("koefficient_of_melting_meteorite")

energy_to_heating_meteorite = thermal_energy_law.law.subs({
    thermal_energy_law.heat_capacity:
    specific_qty_law.law.rhs.subs({
    specific_qty_law.specific_quantity: specific_heat_heating_meteorite,
    specific_qty_law.mass: mass_of_meteorite,
    }),
    thermal_energy_law.temperature_origin:
        temperature_of_meteorite,
    thermal_energy_law.temperature_end:
        temperature_of_meteorite_melting
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
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_before):
    koefficient_of_transfer_to_kinetic * kinetic_energy_of_meteorite,
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_after):
    energy_to_heating_meteorite + energy_to_meteorite_melting
})
mass_of_melting_value = solve(energy_conservation_equation, mass_of_melting_meteorite,
    dict=True)[0][mass_of_melting_meteorite]
koefficient_of_melting_meteorite_value = mass_of_melting_value / mass_of_meteorite
answer = Eq(mass_of_melting_meteorite / mass_of_meteorite, koefficient_of_melting_meteorite_value)
print(f"Total equation:\n{print_expression(answer)}")

koefficient_of_melting_meteorite_percents = koefficient_of_melting_meteorite_value.subs({
    velocity_of_meteorite: Quantity(1.5 * units.kilometers / units.second),
    temperature_of_meteorite: Quantity(300 * units.kelvins),
    koefficient_of_transfer_to_kinetic: Quantity(80 * units.percents),
    specific_heat_heating_meteorite: Quantity(460 * units.joules / (units.kilogram * units.kelvin)),
    specific_heat_melting_meteorite: Quantity(270 * prefixes.kilo * units.joules / units.kilogram),
    temperature_of_meteorite_melting: to_kelvin_quantity(Celsius(1400)),
})
print(
    f"The fraction of a meteorite that melts is: {print_expression(convert_to(Quantity(koefficient_of_melting_meteorite_percents), units.percents).evalf(5))} %"
)
