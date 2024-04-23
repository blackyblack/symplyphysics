#!/usr/bin/env python3
from sympy import solve, symbols
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import units
from symplyphysics.laws.thermodynamics import temperature_is_constant as isothermal_law
from symplyphysics.laws.thermodynamics import zero_heat_transfer as adiabatic_law
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

solved_isothermal_law = solve(isothermal_law.derived_law,
    (isothermal_law.pressure_start, isothermal_law.temperature_end, isothermal_law.pressure_end),
    dict=True)[0][isothermal_law.pressure_end]

solved_adiabatic_law = solve(adiabatic_law.law,
    (adiabatic_law.pressure_start, adiabatic_law.temperature_end, adiabatic_law.pressure_end),
    dict=True)[0][adiabatic_law.pressure_end]

# We need adiabatic volume law to detect the point on the 'volume' axis where
# we adiabatic compression should begin in order to close the cycle
solved_adiabatic_volume_law = solve(adiabatic_law.law,
    (adiabatic_law.pressure_start, adiabatic_law.pressure_end, adiabatic_law.volume_end),
    dict=True)[0][adiabatic_law.volume_end]

carnot_cycle_volume = symbols("carnot_cycle_volume")

GAS_MOLE_COUNT = 1
GAS_TEMPERATURE_START = 200
GAS_TEMPERATURE_END = 100
GAS_SPECIFIC_HEATS_RATIO = 1.66
GAS_VOLUME_START = 1
GAS_VOLUME_ADIABATIC_START = GAS_VOLUME_START + 1

gas_volume_adiabatic_end = solved_adiabatic_volume_law.subs({
    adiabatic_law.temperature_start: GAS_TEMPERATURE_START,
    adiabatic_law.volume_start: GAS_VOLUME_ADIABATIC_START,
    adiabatic_law.specific_heats_ratio: GAS_SPECIFIC_HEATS_RATIO,
    adiabatic_law.temperature_end: GAS_TEMPERATURE_END
})

# Use reversed adiabatic process here: looks like adiabatic expansion
# but we use it to calculate adiabatic compression
gas_volume_isothermal_end = solved_adiabatic_volume_law.subs({
    adiabatic_law.temperature_start: GAS_TEMPERATURE_START,
    adiabatic_law.volume_start: GAS_VOLUME_START,
    adiabatic_law.specific_heats_ratio: GAS_SPECIFIC_HEATS_RATIO,
    adiabatic_law.temperature_end: GAS_TEMPERATURE_END
})

result_pressure_isothermal_expansion = solved_isothermal_law.subs({
    isothermal_law.temperature_start: GAS_TEMPERATURE_START,
    units.molar_gas_constant: 1,
    ideal_gas_equation.mole_count: GAS_MOLE_COUNT,
    isothermal_law.volume_end: carnot_cycle_volume
})

result_pressure_isothermal_compression = solved_isothermal_law.subs({
    isothermal_law.temperature_start: GAS_TEMPERATURE_END,
    units.molar_gas_constant: 1,
    ideal_gas_equation.mole_count: GAS_MOLE_COUNT,
    isothermal_law.volume_end: carnot_cycle_volume
})

result_pressure_adiabatic_expansion = solved_adiabatic_law.subs({
    adiabatic_law.temperature_start: GAS_TEMPERATURE_START,
    adiabatic_law.volume_start: GAS_VOLUME_ADIABATIC_START,
    units.molar_gas_constant: 1,
    ideal_gas_equation.mole_count: GAS_MOLE_COUNT,
    adiabatic_law.specific_heats_ratio: GAS_SPECIFIC_HEATS_RATIO,
    adiabatic_law.volume_end: carnot_cycle_volume
})

result_pressure_adiabatic_compression = solved_adiabatic_law.subs({
    adiabatic_law.temperature_start: GAS_TEMPERATURE_END,
    adiabatic_law.volume_start: gas_volume_isothermal_end,
    units.molar_gas_constant: 1,
    ideal_gas_equation.mole_count: GAS_MOLE_COUNT,
    adiabatic_law.specific_heats_ratio: GAS_SPECIFIC_HEATS_RATIO,
    adiabatic_law.volume_end: carnot_cycle_volume
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
