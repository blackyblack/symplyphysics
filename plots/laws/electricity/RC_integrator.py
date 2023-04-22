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
from symplyphysics.laws.electricity.circuits import resistor_and_capacitor_as_integrator_node as rc_node

time = symbols("time")

initial_voltage = 1
example_capacitance = 1
example_impedance = 1

rc_time_constant = example_capacitance * example_impedance

applied_law = rc_node.law.subs({
    rc_node.time: time,
    rc_node.resistance: example_impedance,
    rc_node.capacitance: example_capacitance,
    rc_node.initial_voltage: initial_voltage
})
capacitor_voltage_function = solve(
    applied_law, rc_node.capacitor_voltage(time),
    dict=True)[0][rc_node.capacitor_voltage(time)]

# see resistor_and_capacitor_as_integrator_node.capacitor_current_eq for a proof the current on resistor equals to current on capacitor
# see resistor_and_capacitor_as_integrator_node.resistor_voltage_eq for a proof the voltage on resistor is (initial_voltage - capacitor_voltage_function)
# see symplyphysics.laws.electricity.current_is_proportional_to_voltage for a proof that current through resistor = resistor voltage / resistor impedance
capacitor_current_function = (initial_voltage -
                              capacitor_voltage_function) / example_impedance

UC = plot(capacitor_voltage_function, (time, 0, 8 * rc_time_constant),
          line_color="blue",
          title="Capacitor voltage",
          label="Capacitor voltage",
          legend=True,
          annotations={},
          backend=MatplotlibBackend,
          show=False)

IC = plot(capacitor_current_function, (time, 0, 8 * rc_time_constant),
          line_color="orange",
          label="Capacitor current",
          legend=True,
          backend=MatplotlibBackend,
          show=False)
UC.append(IC[0])

voltage063 = plot(0.63 * initial_voltage, (time, 0, rc_time_constant),
                  line_color="yellow",
                  label="U = 0.63 of U0",
                  backend=MatplotlibBackend,
                  show=False)
UC.append(voltage063[0])

tau_line = plot(
    1000 * (time - rc_time_constant) * 0.63 * initial_voltage,
    (time, rc_time_constant, rc_time_constant + 0.001),
    label="time = Tau",
    line_color="yellow",
    show=False,
)
UC.append(tau_line[0])

voltage095 = plot(0.95 * initial_voltage, (time, 0, 3 * rc_time_constant),
                  line_color="green",
                  label="U = 0.95 of U0",
                  backend=MatplotlibBackend,
                  show=False)
UC.append(voltage095[0])

tau3_line = plot(
    1000 * (time - 3 * rc_time_constant) * 0.95 * initial_voltage,
    (time, 3 * rc_time_constant, 3 * rc_time_constant + 0.001),
    label="time = 3 * Tau",
    line_color="green",
    show=False,
)
UC.append(tau3_line[0])

voltageFull = plot(initial_voltage, (time, 0, 8 * rc_time_constant),
                   line_color="red",
                   label="U0",
                   backend=MatplotlibBackend,
                   show=False)
UC.append(voltageFull[0])

UC.show()
