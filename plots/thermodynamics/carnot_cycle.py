#!/usr/bin/env python3
from sympy import Eq, solve, symbols
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import units
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.laws.thermodynamics import pressure_and_volume_in_isothermal_process as isothermal_law
from symplyphysics.laws.thermodynamics import pressure_and_volume_in_adiabatic_process as adiabatic_law
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

_temperature_start = symbols("temperature_start")
_temperature_end = symbols("temperature_end")

_isothermal_condition = Eq(_temperature_start, _temperature_end)

ideal_gas_equation_start = ideal_gas_equation.law.subs({
    ideal_gas_equation.pressure: isothermal_law.initial_pressure,
    ideal_gas_equation.volume: isothermal_law.initial_volume,
    ideal_gas_equation.temperature: _temperature_start,
})

ideal_gas_equation_end = ideal_gas_equation.law.subs({
    ideal_gas_equation.pressure: isothermal_law.final_pressure,
    ideal_gas_equation.volume: isothermal_law.final_volume,
    ideal_gas_equation.temperature: _temperature_end,
})

isothermal_pressure = solve(
    [ideal_gas_equation_start, ideal_gas_equation_end, _isothermal_condition],
    (isothermal_law.initial_pressure, _temperature_end, isothermal_law.final_pressure),
    dict=True,
)[0][isothermal_law.final_pressure]
isothermal_pressure = evaluate_expression(isothermal_pressure)

adiabatic_pressure = solve(
    adiabatic_law.law,
    (adiabatic_law.initial_pressure, adiabatic_law.final_temperature, adiabatic_law.final_pressure),
    dict=True,
)[0][adiabatic_law.final_pressure]
adiabatic_pressure = evaluate_expression(adiabatic_pressure)

# We need adiabatic volume law to detect the point on the 'volume' axis where
# we adiabatic compression should begin in order to close the cycle
adiabatic_volume = solve(
    adiabatic_law.law,
    (adiabatic_law.initial_pressure, adiabatic_law.final_pressure, adiabatic_law.final_volume),
    dict=True,
)[0][adiabatic_law.final_volume]
adiabatic_volume = evaluate_expression(adiabatic_volume)

carnot_cycle_volume = symbols("carnot_cycle_volume")

GAS_MOLE_COUNT = 1
GAS_TEMPERATURE_START = 200
GAS_TEMPERATURE_END = 100
GAS_SPECIFIC_HEATS_RATIO = 1.66
GAS_VOLUME_START = 1
GAS_VOLUME_ADIABATIC_START = GAS_VOLUME_START + 1

gas_volume_adiabatic_end = adiabatic_volume.subs({
    adiabatic_law.initial_temperature: GAS_TEMPERATURE_START,
    adiabatic_law.initial_volume: GAS_VOLUME_ADIABATIC_START,
    adiabatic_law.adiabatic_index: GAS_SPECIFIC_HEATS_RATIO,
    adiabatic_law.final_temperature: GAS_TEMPERATURE_END,
    ideal_gas_equation.amount_of_substance: GAS_MOLE_COUNT,
})

# Use reversed adiabatic process here: looks like adiabatic expansion
# but we use it to calculate adiabatic compression
gas_volume_isothermal_end = adiabatic_volume.subs({
    adiabatic_law.initial_temperature: GAS_TEMPERATURE_START,
    adiabatic_law.initial_volume: GAS_VOLUME_START,
    adiabatic_law.adiabatic_index: GAS_SPECIFIC_HEATS_RATIO,
    adiabatic_law.final_temperature: GAS_TEMPERATURE_END,
    ideal_gas_equation.amount_of_substance: GAS_MOLE_COUNT,
})

result_pressure_isothermal_expansion = isothermal_pressure.subs({
    _temperature_start: GAS_TEMPERATURE_START,
    units.molar_gas_constant: 1,
    ideal_gas_equation.amount_of_substance: GAS_MOLE_COUNT,
    isothermal_law.final_volume: carnot_cycle_volume
})

result_pressure_isothermal_compression = isothermal_pressure.subs({
    _temperature_start: GAS_TEMPERATURE_END,
    units.molar_gas_constant: 1,
    ideal_gas_equation.amount_of_substance: GAS_MOLE_COUNT,
    isothermal_law.final_volume: carnot_cycle_volume
})

result_pressure_adiabatic_expansion = adiabatic_pressure.subs({
    adiabatic_law.initial_temperature: GAS_TEMPERATURE_START,
    adiabatic_law.initial_volume: GAS_VOLUME_ADIABATIC_START,
    units.molar_gas_constant: 1,
    ideal_gas_equation.amount_of_substance: GAS_MOLE_COUNT,
    adiabatic_law.adiabatic_index: GAS_SPECIFIC_HEATS_RATIO,
    adiabatic_law.final_volume: carnot_cycle_volume
})

result_pressure_adiabatic_compression = adiabatic_pressure.subs({
    adiabatic_law.initial_temperature: GAS_TEMPERATURE_END,
    adiabatic_law.initial_volume: gas_volume_isothermal_end,
    units.molar_gas_constant: 1,
    ideal_gas_equation.amount_of_substance: GAS_MOLE_COUNT,
    adiabatic_law.adiabatic_index: GAS_SPECIFIC_HEATS_RATIO,
    adiabatic_law.final_volume: carnot_cycle_volume
})

p1 = plot(result_pressure_isothermal_expansion,
    (carnot_cycle_volume, GAS_VOLUME_START, GAS_VOLUME_ADIABATIC_START),
    line_color="blue",
    title="Carnot Cycle",
    xlabel="Volume",
    ylabel="Pressure",
    label="Isothermal",
    legend=True,
    backend=MatplotlibBackend,
    show=False)

p2 = plot(result_pressure_adiabatic_expansion,
    (carnot_cycle_volume, GAS_VOLUME_ADIABATIC_START, gas_volume_adiabatic_end),
    line_color="red",
    label="Adiabatic",
    backend=MatplotlibBackend,
    show=False)

p3 = plot(result_pressure_isothermal_compression,
    (carnot_cycle_volume, gas_volume_isothermal_end, gas_volume_adiabatic_end),
    line_color="blue",
    backend=MatplotlibBackend,
    label="",
    show=False)

p4 = plot(result_pressure_adiabatic_compression,
    (carnot_cycle_volume, GAS_VOLUME_START, gas_volume_isothermal_end),
    line_color="red",
    backend=MatplotlibBackend,
    label="",
    show=False)

p1.append(p2[0])
p1.append(p3[0])
p1.append(p4[0])
p1[2].label = ""
p1[3].label = ""

p1.show()
