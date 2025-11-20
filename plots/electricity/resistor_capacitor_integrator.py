#!/usr/bin/env python3
r"""
For an :math:`RC` integrator (i.e. a serial :math:`RC` circuit), assume a step voltage input and
plot the voltage across charging capacitor and the current in the circuit as functions of time.

The step input is defined as
.. math::
    U(t) = \begin{cases} 0 & t < 0 \\ U_0 & t \\ge 0 \end{cases}

Note that the time constant of the serial :math:`RC` circuit is :math:`\tau = R C`. After time
:math:`t = N \tau` passes, the voltage across the capacitor will reach :math:`1 - e^{-N}` times
its maximum value. The circuit current, on the other hand, will simply decrease by a factor of
:math:`e^N`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/RC_circuit#Time-domain_considerations>`__.
"""

from sympy import symbols as sym_symbols, solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import units, convert_to_si
from symplyphysics.laws.electricity.circuits import (
    voltage_across_charging_capacitor_in_serial_resistor_capacitor_circuit as rc_node,
    time_constant_of_resistor_capacitor_circuit as time_constant_law,
)

time = sym_symbols("time")

INITIAL_VOLTAGE = convert_to_si(1 * units.volt)
EXAMPLE_CAPACITANCE = convert_to_si(1 * units.farad)
EXAMPLE_RESISTANCE = convert_to_si(1 * units.ohm)

RC_TIME_CONSTANT = EXAMPLE_CAPACITANCE * EXAMPLE_RESISTANCE

time_constant_expr = time_constant_law.law.rhs.subs({
    time_constant_law.resistance: EXAMPLE_RESISTANCE,
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
capacitor_current_function = (INITIAL_VOLTAGE - capacitor_voltage_function) / EXAMPLE_RESISTANCE

UC = plot(capacitor_voltage_function, (time, 0, 8 * RC_TIME_CONSTANT),
    line_color="blue",
    title="Charging capacitor in a serial $RC$ circuit",
    label="Capacitor voltage",
    xlabel="time $t$",
    ylabel="reduced quantities",
    legend=True,
    backend=MatplotlibBackend,
    show=False)

IC = plot(capacitor_current_function, (time, 0, 8 * RC_TIME_CONSTANT),
    line_color="orange",
    label="Circuit current",
    legend=True,
    backend=MatplotlibBackend,
    show=False)
UC.append(IC[0])

voltage063 = plot(0.63 * INITIAL_VOLTAGE, (time, 0, RC_TIME_CONSTANT),
    line_color="yellow",
    label="$U = 0.63 U_0$",
    backend=MatplotlibBackend,
    show=False)
UC.append(voltage063[0])

tau_line = plot(
    1000 * (time - RC_TIME_CONSTANT) * 0.63 * INITIAL_VOLTAGE,
    (time, RC_TIME_CONSTANT, RC_TIME_CONSTANT + 0.001),
    label="$t = \\tau$",
    line_color="yellow",
    show=False,
)
UC.append(tau_line[0])

voltage095 = plot(0.95 * INITIAL_VOLTAGE, (time, 0, 3 * RC_TIME_CONSTANT),
    line_color="green",
    label="$U = 0.95 U_0$",
    backend=MatplotlibBackend,
    show=False)
UC.append(voltage095[0])

tau3_line = plot(
    1000 * (time - 3 * RC_TIME_CONSTANT) * 0.95 * INITIAL_VOLTAGE,
    (time, 3 * RC_TIME_CONSTANT, 3 * RC_TIME_CONSTANT + 0.001),
    label="$t = 3 \\tau$",
    line_color="green",
    show=False,
)
UC.append(tau3_line[0])

voltageFull = plot(INITIAL_VOLTAGE, (time, 0, 8 * RC_TIME_CONSTANT),
    line_color="red",
    label="$U_0$",
    backend=MatplotlibBackend,
    show=False)
UC.append(voltageFull[0])

UC.show()
