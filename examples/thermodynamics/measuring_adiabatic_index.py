#!/usr/bin/env python3

from sympy import solve, Eq, symbols
from symplyphysics import print_expression
from symplyphysics.definitions import heat_capacity_ratio
from symplyphysics.laws.electricity import energy_via_constant_power_and_time as energy_law
from symplyphysics.laws.thermodynamics import heat_is_heat_capacity_times_temperature_change as thermal_eqn
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation
from symplyphysics.laws.quantities import quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law

# Description
## In order to measure the value of the adiabatic index (gamma = C_p / C_V), one can use the
## following method. A certain amount of gas of known initial volume and pressure is heated twice
## with an electric coil for the same amount of time, firstly at constant volume measuring the
## resulting pressure, and secondly at constant pressure measuring the resulting volume. Calculate
## the adiabatic index from this data. Assume that the gas is ideal and that the gas chamber is adiabatically
## isolated from the environment.

amount_of_gas = symbols("amount_of_gas")
initial_pressure, final_isochoric_pressure = symbols("initial_pressure final_isochoric_pressure")
initial_volume, final_isobaric_volume = symbols("initial_volume final_isobaric_volume")

# Calculate temperature values from the ideal gas equation of state

temperature_expr = solve(ideal_gas_equation.law,
    ideal_gas_equation.temperature)[0].subs(ideal_gas_equation.amount_of_substance, amount_of_gas)

initial_temperature_expr = temperature_expr.subs({
    ideal_gas_equation.pressure: initial_pressure,
    ideal_gas_equation.volume: initial_volume,
})

final_isochoric_temperature_expr = temperature_expr.subs({
    ideal_gas_equation.pressure: final_isochoric_pressure,
    ideal_gas_equation.volume: initial_volume,
})

final_isobaric_temperature_expr = temperature_expr.subs({
    ideal_gas_equation.pressure: initial_pressure,
    ideal_gas_equation.volume: final_isobaric_volume,
})

# Calculate amounts of heat in both cases using the temperature values above

isochoric_molar_heat, isobaric_molar_heat = symbols("isochoric_molar_heat isobaric_molar_heat")

isochoric_heat_capacity = molar_qty_law.law.rhs.subs({
    molar_qty_law.molar_quantity: isochoric_molar_heat,
    molar_qty_law.amount_of_substance: amount_of_gas,
})

isochoric_amount_of_heat = thermal_eqn.law.rhs.subs({
    thermal_eqn.heat_capacity: isochoric_heat_capacity,
    thermal_eqn.temperature_change: final_isobaric_temperature_expr - initial_temperature_expr,
}).simplify()

isobaric_heat_capacity = molar_qty_law.law.rhs.subs({
    molar_qty_law.molar_quantity: isobaric_molar_heat,
    molar_qty_law.amount_of_substance: amount_of_gas,
})

isobaric_amount_of_heat = thermal_eqn.law.rhs.subs({
    thermal_eqn.heat_capacity: isobaric_heat_capacity,
    thermal_eqn.temperature_change: final_isobaric_temperature_expr - initial_temperature_expr,
}).simplify()

# Calculate the adiabatic index using the amounts of heat during two heating instances.

solved = solve(
    [
        energy_law.law.subs(energy_law.energy, isochoric_amount_of_heat),
        energy_law.law.subs(energy_law.energy, isobaric_amount_of_heat)
    ],
    (isochoric_molar_heat, isobaric_molar_heat),
    dict=True,
)[0]

isochoric_molar_heat_expr = solved[isochoric_molar_heat]
isobaric_molar_heat_expr = solved[isobaric_molar_heat]

adiabatic_index_expr = heat_capacity_ratio.definition.rhs.subs({
    heat_capacity_ratio.isobaric_heat_capacity: isobaric_molar_heat_expr,
    heat_capacity_ratio.isochoric_heat_capacity: isochoric_molar_heat_expr,
})

# We can show that the adiabatic index depends only on the factors by which
# pressure and volume increase in the experiment

isochoric_pressure_increase_factor = symbols("isochoric_pressure_increase_factor")
isochoric_factor_eqn = Eq(
    final_isochoric_pressure,
    initial_pressure * isochoric_pressure_increase_factor,
)

isobaric_volume_increase_factor = symbols("isobaric_volume_increase_factor")
isobaric_factor_eqn = Eq(
    final_isobaric_volume,
    initial_volume * isobaric_volume_increase_factor,
)

adiabatic_index = symbols("adiabatic_index")
adiabatic_index_expr_ = solve(
    [Eq(adiabatic_index, adiabatic_index_expr), isochoric_factor_eqn, isobaric_factor_eqn],
    (adiabatic_index, initial_pressure, initial_volume),
    dict=True,
)[0][adiabatic_index]

# Display results

print(
    "The adiabatic index of an ideal gas:",
    print_expression(Eq(adiabatic_index, adiabatic_index_expr)),
    "An alternative expression via increase factors:",
    print_expression(Eq(adiabatic_index, adiabatic_index_expr_)),
    sep="\n\n",
)
