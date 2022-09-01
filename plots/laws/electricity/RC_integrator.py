#!/usr/bin/env python3

# Description
## Time constant of RC Tau = R*C (seconds).
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
# from symplyphysics.laws.electricity import sums_of_all_voltages_in_loop_is_zero as kirchhof_law_2

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
example_impedance = 3000

tau = example_capacitance * example_impedance

# capacitor_voltage = capacitance_definition.calculate_voltage()

# U(t) shoud be U0(1 - e**(-t/RC))
# capacitor_voltage = initial_voltage * (1 - exp(-time/(example_capacitance * example_impedance)))

UC_func = initial_voltage * (1 - exp(-time/(example_capacitance * example_impedance)))
IC_func   = (initial_voltage / example_impedance) * exp(-time/(example_capacitance * example_impedance)) * 2500 #2500 is to scale curve of current to voltage plot

UC = plot(
    UC_func,
    (time, 0, 8 * tau),
    line_color='blue',    
    title='Capacitor voltage',    
    label = 'Capacitor voltage',    
    legend=True,    
    annotations = {},
    backend=MatplotlibBackend,    
    show=False)        

IC = plot(
    IC_func,
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
    (time, 0, 6 * tau),
    line_color='red',      
    label='U0',
    backend=MatplotlibBackend,
    show=False)
UC.append(voltageFull[0])

UC.show()
