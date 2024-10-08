#!/usr/bin/env python3

# Description
## Time constant of RC Tau = R*C (seconds).
## RC integrator is a circuit with capacitor and resistor in series. Initial_voltage is applied to whole circuit and integrated voltage is obtained from capacitor.
## If some voltage is applied to RC integrator, capacitor starts to charge and it's voltage rises.
## Capacitor voltage will never reach initial voltage.
## Voltage on the capacitor reaches 63% of initial voltage after 1*Tau seconds, and reaches 95% after 3*Tau seconds
## independently of initial voltage.

from sympy import symbols, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics.laws.electricity.circuits import (
    voltage_across_charging_capacitor_in_serial_resistor_capacitor_circuit as rc_node,
    time_constant_of_resistor_capacitor_circuit as time_constant_law,
)

time = symbols("time")

INITIAL_VOLTAGE = 1
EXAMPLE_CAPACITANCE = 1
EXAMPLE_IMPEDANCE = 1

RC_TIME_CONSTANT = EXAMPLE_CAPACITANCE * EXAMPLE_IMPEDANCE

time_constant_expr = time_constant_law.law.rhs.subs({
    time_constant_law.resistance: EXAMPLE_IMPEDANCE,
    time_constant_law.capacitance: EXAMPLE_CAPACITANCE,
})

applied_law = rc_node.law.subs({
    rc_node.time: time,
    rc_node.time_constant: time_constant_expr,
    rc_node.source_voltage: INITIAL_VOLTAGE
})
capacitor_voltage_function = solve(applied_law, rc_node.capacitor_voltage,
    dict=True)[0][rc_node.capacitor_voltage]

# see resistor_and_capacitor_as_integrator_node.capacitor_current_eq for a proof the current on resistor equals to current on capacitor
# see resistor_and_capacitor_as_integrator_node.resistor_voltage_eq for a proof the voltage on resistor is (initial_voltage - capacitor_voltage_function)
# see symplyphysics.laws.electricity.current_is_proportional_to_voltage for a proof that current through resistor = resistor voltage / resistor impedance
capacitor_current_function = (INITIAL_VOLTAGE - capacitor_voltage_function) / EXAMPLE_IMPEDANCE

UC = plot(capacitor_voltage_function, (time, 0, 8 * RC_TIME_CONSTANT),
    line_color="blue",
    title="Capacitor voltage",
    label="Capacitor voltage",
    legend=True,
    backend=MatplotlibBackend,
    show=False)

IC = plot(capacitor_current_function, (time, 0, 8 * RC_TIME_CONSTANT),
    line_color="orange",
    label="Capacitor current",
    legend=True,
    backend=MatplotlibBackend,
    show=False)
UC.append(IC[0])

voltage063 = plot(0.63 * INITIAL_VOLTAGE, (time, 0, RC_TIME_CONSTANT),
    line_color="yellow",
    label="U = 0.63 of U0",
    backend=MatplotlibBackend,
    show=False)
UC.append(voltage063[0])

tau_line = plot(
    1000 * (time - RC_TIME_CONSTANT) * 0.63 * INITIAL_VOLTAGE,
    (time, RC_TIME_CONSTANT, RC_TIME_CONSTANT + 0.001),
    label="time = Tau",
    line_color="yellow",
    show=False,
)
UC.append(tau_line[0])

voltage095 = plot(0.95 * INITIAL_VOLTAGE, (time, 0, 3 * RC_TIME_CONSTANT),
    line_color="green",
    label="U = 0.95 of U0",
    backend=MatplotlibBackend,
    show=False)
UC.append(voltage095[0])

tau3_line = plot(
    1000 * (time - 3 * RC_TIME_CONSTANT) * 0.95 * INITIAL_VOLTAGE,
    (time, 3 * RC_TIME_CONSTANT, 3 * RC_TIME_CONSTANT + 0.001),
    label="time = 3 * Tau",
    line_color="green",
    show=False,
)
UC.append(tau3_line[0])

voltageFull = plot(INITIAL_VOLTAGE, (time, 0, 8 * RC_TIME_CONSTANT),
    line_color="red",
    label="U0",
    backend=MatplotlibBackend,
    show=False)
UC.append(voltageFull[0])

UC.show()
