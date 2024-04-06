#!/usr/bin/env python3

from sympy import solve, Eq, dsolve
from symplyphysics import print_expression, Quantity, units, convert_to, Symbol
from symplyphysics.definitions import power_is_energy_derivative as power_def
from symplyphysics.laws.thermodynamics import (
    thermal_energy_from_mass_and_temperature as heat_law,
    pressure_from_temperature_and_volume as ideal_gas_law,
)

# Description
## In order to measure the value of the adiabatic index (gamma = C_p / C_V), one can use the following method. A certain amount of gas of known initial volume and pressure is heated twice with an electric coil for the same amount of time, firstly at constant volume measuring the resulting pressure, and secondly at constant pressure measuring the resulting volume. Calculate the adiabatic index from this data. Assume that the gas is ideal.

pprint = lambda *args, **kwargs: print(*map(print_expression, args), **kwargs)

gas_mass = Symbol("gas_mass", units.mass)
initial_pressure = Symbol("initial_pressure", units.pressure)
initial_volume = Symbol("initial_volume", units.volume)
final_pressure = Symbol("final_pressure", units.pressure)
final_volume = Symbol("final_volume", units.volume)

isochoric_specific_heat_capacity = Symbol(
    "isochoric_specific_heat_capacity",
    units.energy / (units.mass * units.temperature),
)
isobaric_specific_heat_capacity = Symbol(
    "isobaric_specific_heat_capacity",
    units.energy / (units.mass * units.temperature),
)

# The power of the coil is constant during both instances of heating.
heating_power = Symbol("heating_power", units.power)

power_eqn = power_def.definition.subs({
    power_def.power(power_def.time): heating_power,
})

heat_eqn = dsolve(
    power_eqn,
    power_def.energy(power_def.time),
    ics={power_def.energy(0): 0},  # heat transferred is zero when no heating has started yet
)

# TODO: write about adiabatic isolation

temperature_difference = Symbol("temperature_difference", units.temperature)

isochoric_heat = heat_law.law.rhs.subs({
    heat_law.specific_heat_capacity: isochoric_specific_heat_capacity,
    heat_law.temperature_end: temperature_difference,
    heat_law.temperature_origin: 0,
})

isobaric_heat = heat_law.law.rhs.subs({
    heat_law.specific_heat_capacity: isobaric_specific_heat_capacity,
    heat_law.temperature_end: temperature_difference,
    heat_law.temperature_origin: 0,
})

# Since the time of heating is the same, we can write the following set of equations

equal_amounts_of_heat = [
    heat_eqn.subs(power_def.energy(power_def.time), isochoric_heat),
    heat_eqn.subs(power_def.energy(power_def.time), isobaric_heat),
]


