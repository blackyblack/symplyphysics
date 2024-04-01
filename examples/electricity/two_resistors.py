#!/usr/bin/env python3

from sympy import Idx, solve, symbols
from symplyphysics import units, global_index
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohm_law
from symplyphysics.laws.electricity.circuits import conductivity_of_two_parallel_resistors as parallel_conductance
from symplyphysics.laws.electricity.circuits import resistivity_of_serial_resistors as serial_resistance
from symplyphysics.definitions import electrical_conductivity_is_inversed_resistance as conductivity_law

# Two resistors are connected across a 12 V battery with internal resistance of 1 Ohm.
# When they are connected in parallel, the current in the circuit is 4 A.
# When they are connected in series, the current in the circuit is 1 A.
# Find the resistance of said resistors.

R1 = symbols("R1")
R2 = symbols("R2")
R_battery = symbols("R_battery")
E_battery = symbols("E_battery")
I_series = symbols("I_series")
I_parallel = symbols("I_parallel")

# Parallel connection

## Find resistance R12 using the law of conductance

sigma1 = solve(conductivity_law.definition,
    conductivity_law.object_conductivity)[0].subs({conductivity_law.object_resistance: R1})

sigma2 = solve(conductivity_law.definition,
    conductivity_law.object_conductivity)[0].subs({conductivity_law.object_resistance: R2})

parallel_law = parallel_conductance.law.subs({
    parallel_conductance.first_conductance: sigma1,
    parallel_conductance.second_conductance: sigma2,
})
sigma_parallel = solve(parallel_law, parallel_conductance.parallel_conductance)[0]

resistance_definition = conductivity_law.definition.subs(
    {conductivity_law.object_conductivity: sigma_parallel})
R12_parallel = solve(resistance_definition, conductivity_law.object_resistance)[0]

## The sum resistance R12 is still connected in series to the internal resistance of the battery

index_local = Idx("index_local", (1, 2))
parallel_law_two_resistors = serial_resistance.law.subs(global_index, index_local).doit()

total_parallel_law = parallel_law_two_resistors.subs({
    serial_resistance.resistance[1]: R_battery,
    serial_resistance.resistance[2]: R12_parallel,
})
R_total_parallel = solve(total_parallel_law, serial_resistance.serial_resistance)[0]

ohm_law_parallel = ohm_law.law.subs({
    ohm_law.current: I_parallel,
    ohm_law.voltage: E_battery,
    ohm_law.resistance: R_total_parallel,
})

# Serial connection

serial_law_two_resistors = serial_resistance.law.subs(global_index, index_local).doit()

serial_law = serial_law_two_resistors.subs({
    serial_resistance.resistance[1]: R1,
    serial_resistance.resistance[2]: R2,
})
R12_serial = solve(serial_law, serial_resistance.serial_resistance)[0]

serial_law_resistors_and_battery = serial_resistance.law.subs(global_index, index_local).doit()

total_serial_law = serial_law_resistors_and_battery.subs({
    serial_resistance.resistance[1]: R_battery,
    serial_resistance.resistance[2]: R12_serial,
})
R_total_serial = solve(total_serial_law, serial_resistance.serial_resistance)[0]

ohm_law_serial = ohm_law.law.subs({
    ohm_law.current: I_series,
    ohm_law.voltage: E_battery,
    ohm_law.resistance: R_total_serial,
})

# Solve the system of two equations resulting from the application of Ohm's law
# There is also a second solution that gives the same results in a reverse order
result_R1, result_R2 = solve([ohm_law_parallel, ohm_law_serial], R1, R2)[0]

known_values = {R_battery: 1, E_battery: 12, I_series: 1, I_parallel: 4}

result_R1_value = result_R1.subs(known_values).evalf(3)
result_R2_value = result_R2.subs(known_values).evalf(3)

# The resulting formulas are too long for printing
print(f"R1 = {result_R1_value} {units.ohm}, R2 = {result_R2_value} {units.ohm}")
