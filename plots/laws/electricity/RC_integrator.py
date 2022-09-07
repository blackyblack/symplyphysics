#!/usr/bin/env python3

# Description
## Time constant of RC Tau = R*C (seconds).
## RC integrator is a circuit with capacitor and resistor in series. Initial_voltage is applied to whole circuit and integrated voltage is obtained from capacitor.
## If some voltage is applied to RC integrator, capacitor starts to charge and it's voltage rises. 
## Capacitor voltage will never reach initial voltage.
## Voltage on the capacitor reaches 63% of initial voltage after 1*Tau seconds, and reaches 95% after 3*Tau seconds
## independently of initial voltage.

from sympy import solve, dsolve, symbols, Function, Eq, exp
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics.definitions import current_is_charge_derivative as charge_definition
from symplyphysics.definitions import capacitance_from_charge_and_voltage as capacitance_definition
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohms_law

# TODO: prove that ohms_law.voltage = voltage_initial - voltage_function(time)
# According to Kirchhoff's laws #1 and #2 we have:
# 1. capacitor_current(time) = resistor_current(time)
# 2. initial_voltage = capacitor_voltage(time) + resistor_voltage_time
# According to Ohm's law we have:
# 3. resistor_current(time) = resistor_voltage(time) / resistor_impedance (R)
# According to capacitance definition we have:
# 4. capacitor_charge(time) = capacitor_voltage(time) * capacitor_capacitance (C)
# According to current definition:
# 5. capacitor_current(time) = d(capacitor_charge)/d(time) OR capacitor_charge(time) = Integral(capacitor_current)d(time)

time = symbols('time')
initial_voltage, resistance = symbols('voltage_initial resistance')
capacitor_voltage_function = symbols('capacitor_voltage_function', cls = Function)

capacitance_eq = capacitance_definition.definition.subs({capacitance_definition.charge: charge_definition.charge_function(time), capacitance_definition.voltage: capacitor_voltage_function(time)})
charge_eq = charge_definition.definition.subs(charge_definition.time, time)

current_eq = ohms_law.law.subs({ohms_law.voltage: (initial_voltage - capacitor_voltage_function(time)), ohms_law.resistance: resistance, ohms_law.current: charge_definition.current_function(time)})

derived_law = [capacitance_eq, charge_eq, current_eq]
solved_charge_function = solve(derived_law, (charge_definition.current_function(time), capacitor_voltage_function(time), charge_definition.charge_function(time)), dict=True)[0][charge_definition.charge_function(time)]

charge_diff_eq = Eq(charge_definition.charge_function(time), solved_charge_function)

C1 = symbols('C1')
charge_on_capacitor_eq = dsolve(charge_diff_eq).subs(C1, 1)
charge_on_capacitor_function = solve(charge_on_capacitor_eq, charge_definition.charge_function(time), dict=True)[0][charge_definition.charge_function(time)]

initial_voltage = 1
example_capacitance = 1
example_impedance = 1

tau = example_capacitance * example_impedance

capacitor_voltage_function_example = initial_voltage * (1 - exp(-time/(example_capacitance * example_impedance)))
capacitor_current_function_example = (initial_voltage - capacitor_voltage_function_example) / example_impedance

UC = plot(
    capacitor_voltage_function_example,
    (time, 0, 8 * tau),
    line_color='blue',    
    title='Capacitor voltage',    
    label = 'Capacitor voltage',    
    legend=True,    
    annotations = {},
    backend=MatplotlibBackend,    
    show=False)        

IC = plot(
    capacitor_current_function_example,
    (time, 0, 8 * tau),
    line_color='orange',    
    label='Capacitor current',    
    legend=True,    
    backend=MatplotlibBackend,
    show=False)
UC.append(IC[0])    

voltage063 = plot(
    0.63 * initial_voltage,
    (time, 0, tau),
    line_color='yellow',        
    label='U = 0.63 of U0',    
    backend=MatplotlibBackend,
    show=False)
UC.append(voltage063[0])    

tau_line = plot(
    1000 * (time - tau) * 0.63 * initial_voltage,
    (time, tau, tau + 0.001),
    label = 'time = tau',
    line_color='yellow',
    show=False,
    )
UC.append(tau_line[0])

voltage095 = plot(
    0.95 * initial_voltage,
    (time, 0, 3 * tau),
    line_color='green',     
    label='U = 0.95 of U0',
    backend=MatplotlibBackend,
    show=False)
UC.append(voltage095[0])

tau3_line = plot(
    1000 * (time - 3 * tau) * 0.95 * initial_voltage,
    (time, 3 * tau, 3 * tau + 0.001),
    label = 'time = 3 * tau',
    line_color='green',
    show=False,
    )
UC.append(tau3_line[0])

voltageFull = plot(
    initial_voltage,
    (time, 0, 8 * tau),
    line_color='red',      
    label='U0',
    backend=MatplotlibBackend,
    show=False)
UC.append(voltageFull[0])

UC.show()
