#!/usr/bin/env python3

from sympy import solve
from symplyphysics import Symbol, units
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohm_law
from symplyphysics.laws.electricity.circuits import resistivity_of_parallel_resistors as parallel_resistors
from symplyphysics.laws.electricity.circuits import resistivity_of_serial_resistors as serial_resistors

# Two resistors are connected across a 12 V battery with internal resistance of 1 Ohm.
# When they are connected in parallel, the current in the circuit is 4 A.
# When they are connected in series, the current in the circuit is 1 A.
# Find the resistance of said resistors.

R1 = Symbol("R1", units.impedance)
R2 = Symbol("R2", units.impedance)
R_battery = Symbol("Rb", units.impedance)
E_battery = Symbol("E", units.voltage)
I_series = Symbol("Is", units.current)
I_parallel = Symbol("Ip", units.current)

# Parallel connection

parallel_law = parallel_resistors.law.subs({parallel_resistors.inv_resistances: (1 / R1, 1 / R2)})
R12_parallel = solve(parallel_law, parallel_resistors.parallel_resistance)[0]

# The sum resistance R12 is still connected in series to the internal resistance of the battery
total_parallel_law = serial_resistors.law.subs({serial_resistors.resistances: (R_battery, R12_parallel)})
R_total_parallel = solve(total_parallel_law, serial_resistors.serial_resistance)[0]

ohm_law_parallel = ohm_law.law.subs({
    ohm_law.current: I_parallel,
    ohm_law.voltage: E_battery,
    ohm_law.resistance: R_total_parallel,
})

# Serial connection

serial_law = serial_resistors.law.subs(serial_resistors.resistances, (R1, R2))
R12_serial = solve(serial_law, serial_resistors.serial_resistance)[0]

total_serial_law = serial_resistors.law.subs(serial_resistors.resistances, (R_battery, R12_serial))
R_total_serial = solve(total_serial_law, serial_resistors.serial_resistance)[0]

ohm_law_serial = ohm_law.law.subs({
    ohm_law.current: I_series,
    ohm_law.voltage: E_battery,
    ohm_law.resistance: R_total_serial,
})

# Solve the system of two equations resulting from the application of Ohm's law
# There is also a second solution that gives the same results in a reverse order
result_R1, result_R2 = solve([ohm_law_parallel, ohm_law_serial], R1, R2)[0]

known_values = {
    R_battery: 1,
    E_battery: 12,
    I_series: 1,
    I_parallel: 4
}

result_R1_value = result_R1.subs(known_values).evalf(3)
result_R2_value = result_R2.subs(known_values).evalf(3)

# The resulting formulas are too long for printing
print(f"R1 = {result_R1_value} {units.ohm}, R2 = {result_R2_value} {units.ohm}")
