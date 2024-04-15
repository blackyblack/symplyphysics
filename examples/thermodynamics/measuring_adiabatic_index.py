#!/usr/bin/env python3

from sympy import solve, Eq, dsolve, symbols
from symplyphysics import print_expression, Quantity, units, convert_to
from symplyphysics.definitions import power_is_energy_derivative as power_def
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as heat_law
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

# Description
## In order to measure the value of the adiabatic index (gamma = C_p / C_V), one can use the following method. A certain amount of gas of known initial volume and pressure is heated twice with an electric coil for the same amount of time, firstly at constant volume measuring the resulting pressure, and secondly at constant pressure measuring the resulting volume. Calculate the adiabatic index from this data. Assume that the gas is ideal.

pprint = lambda *args, **kwargs: print(*map(print_expression, args), **kwargs)

gas_amount = symbols("gas_amount")  # units.amount_of_substance

initial_pressure, final_pressure = symbols("initial_pressure final_pressure")
initial_volume = symbols("initial_volume final_volume")

isochoric_specific_heat = symbols("isochoric_specific_heat")
isobaric_specific_heat = symbols("isobaric_specific_heat")

# The power of the coil is constant during both instances of heating.
heating_power = symbols("heating_power")

power_eqn = power_def.definition.subs({
    power_def.power(power_def.time): heating_power,
})

heat_eqn = dsolve(
    power_eqn,
    power_def.energy(power_def.time),
    ics={power_def.energy(0): 0},  # heat transferred is zero when no heating has started yet
)

# TODO: write about adiabatic isolation

temperature_difference = symbols("temperature_difference")

isochoric_heat = heat_law.law.rhs.subs({
    heat_law.specific_heat_capacity: isochoric_specific_heat,
    heat_law.temperature_end: temperature_difference,
    heat_law.temperature_origin: 0,
})

isobaric_heat = heat_law.law.rhs.subs({
    heat_law.specific_heat_capacity: isobaric_specific_heat,
    heat_law.temperature_end: temperature_difference,
    heat_law.temperature_origin: 0,
})

# Since the time of heating is the same, we can write the following set of equations

equal_amounts_of_heat = [
    heat_eqn.subs(power_def.energy(power_def.time), isochoric_heat),
    heat_eqn.subs(power_def.energy(power_def.time), isobaric_heat),
]


