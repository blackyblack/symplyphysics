#!/usr/bin/env python3

# Description
## Time constant of RC Tau = R*C (seconds).
## If some voltage is applied to RC integrator, capacitor starts to charge and it's voltage rises. 
## Voltage on the capacitor reaches 63% of initial voltage after 1*Tau seconds, and reaches 95% after 3*Tau seconds
## independently of initial voltage.
## For example, MCU is power-supplied with 3.3V. It's Reset pin is pulled up to power rail with internal resistor of 48kOhm
## and externally bypassed to ground with 0.1uF capacitor. MCU starts after Reset pin reaches 2.5V.
## We'll try to estimate this period of time with lelp of diagrams.

from sympy import solve, dsolve, symbols, Function, Eq
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics.definitions import current_is_charge_derivative as charge_definition
from symplyphysics.definitions import capacitance_from_charge_and_voltage as capacitance_definition
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohms_law

time = symbols('time')
voltage_initial, resistance = symbols('voltage_initial resistance')

# derive charge function for the RC integrator node

voltage_function = symbols('voltage_function', cls = Function)

capacitance_eq = capacitance_definition.definition.subs({capacitance_definition.charge: charge_definition.charge_function(time), capacitance_definition.voltage: voltage_function(time)})
charge_eq = charge_definition.definition.subs(charge_definition.time, time)

# TODO: prove that ohms_law.voltage = voltage_initial - voltage_function(time)
current_eq = ohms_law.law.subs({ohms_law.voltage: (voltage_initial - voltage_function(time)), ohms_law.resistance: resistance, ohms_law.current: charge_definition.current_function(time)})

derived_law = [capacitance_eq, charge_eq, current_eq]
solved_charge_function = solve(derived_law, (charge_definition.current_function(time), voltage_function(time), charge_definition.charge_function(time)), dict=True)[0][charge_definition.charge_function(time)]

charge_diff_eq = Eq(charge_definition.charge_function(time), solved_charge_function)

C1 = symbols('C1')
charge_on_capacitor_eq = dsolve(charge_diff_eq).subs(C1, 1)
charge_on_capacitor_function = solve(charge_on_capacitor_eq, charge_definition.charge_function(time), dict=True)[0][charge_definition.charge_function(time)]

initial_voltage = 3.3
example_capacitance = 0.0001
example_impedance = 48000

tau = example_capacitance * example_impedance



capacitor_voltage = capacitance_definition.calculate_voltage()

tau_lim1 = symbols('tau_lim1')
tau_lim3 = symbols('tau_lim3')

p1 = plot(
    capacitor_voltage(time),
    (time, 0, 10),
    line_color='blue',
    title='Capacitor voltage',    
    legend=True,
    backend=MatplotlibBackend,
    show=False)

'''
voltage063 = plot(
    0.63 * initial_voltage,
    (time, 0, 10 * tau),
    line_color='red',    
    backend=MatplotlibBackend,
    show=False)
p1.append(voltage063[0])

voltage095 = plot(
    0.95 * initial_voltage,
    (time, 0, 10 * tau),
    line_color='red',    
    backend=MatplotlibBackend,
    show=False)
p1.append(voltage095[0])
'''
p1.show()
