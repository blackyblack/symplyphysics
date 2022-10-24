#!/usr/bin/env python3

# Description
## You throw the grenade to enemies. Which angle of throwing should you choose to hit the farest enemy?
## Flying velocity V0 vector is the sum of two orthogonal - vertical and horizontal speed vectors V_vert and V_hor.
## If alfa is throwing angle, V_vert = V0 * sin(alpha) and V_hor = V0 * cos(alpha)
## Vertical speed is accelerated by the gravitational acceleration.
## Horizontal speed is constant.

from sympy import sin, cos
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import symbols, Eq, pretty, solve, Function

flight_time = symbols('flight_time')
V0 = symbols('V0')
alpha = symbols('alpha')
Distance = symbols('Distance')
V_vertical = symbols('v_vertical')
V_horizontal = symbols('v_horizontal')
V_vertical = V0 * sin(alpha)
V_horizontal = V0 * cos(alpha)
time = symbols('time')
Length = symbols('Length')
g = 9.806

## Height function of time is h = V_vertical * t - g * t**2 /2
## Flight ends when h == 0

end_of_flight = Eq(0, V_vertical * time - g * time**2 /2)
flight_time = solve(end_of_flight, time, dict=True)[0][time]

print(flight_time) #should be 2 solvations, we need one is not zero
##kostil
flight_time = 2 * V_vertical / g

Length = symbols('Length', cls = Function)
Length = V_horizontal * flight_time(alpha)

print(Length)

Len = plot(
    Length(alpha),
    (alpha, 0, 90),
    line_color='blue',    
    title='Length',    
    label = 'Length',    
    legend=True,    
    annotations = {},
    backend=MatplotlibBackend,    
    show=False)
Len.show()

'''
applied_law = rc_node.law.subs({rc_node.time: time, rc_node.resistance: example_impedance, rc_node.capacitance: example_capacitance, rc_node.initial_voltage: initial_voltage})
capacitor_voltage_function = solve(applied_law, rc_node.capacitor_voltage_function(time), dict=True)[0][rc_node.capacitor_voltage_function(time)]

capacitor_current_function = (initial_voltage - capacitor_voltage_function) / example_impedance

UC = plot(
    capacitor_voltage_function,
    (time, 0, 8 * rc_time_constant),
    line_color='blue',    
    title='Capacitor voltage',    
    label = 'Capacitor voltage',    
    legend=True,    
    annotations = {},
    backend=MatplotlibBackend,    
    show=False)        

IC = plot(
    capacitor_current_function,
    (time, 0, 8 * rc_time_constant),
    line_color='orange',    
    label='Capacitor current',    
    legend=True,    
    backend=MatplotlibBackend,
    show=False)
UC.append(IC[0])    

voltage063 = plot(
    0.63 * initial_voltage,
    (time, 0, rc_time_constant),
    line_color='yellow',        
    label='U = 0.63 of U0',    
    backend=MatplotlibBackend,
    show=False)
UC.append(voltage063[0])    

tau_line = plot(
    1000 * (time - rc_time_constant) * 0.63 * initial_voltage,
    (time, rc_time_constant, rc_time_constant + 0.001),
    label = 'time = Tau',
    line_color='yellow',
    show=False,
    )
UC.append(tau_line[0])

voltage095 = plot(
    0.95 * initial_voltage,
    (time, 0, 3 * rc_time_constant),
    line_color='green',     
    label='U = 0.95 of U0',
    backend=MatplotlibBackend,
    show=False)
UC.append(voltage095[0])

tau3_line = plot(
    1000 * (time - 3 * rc_time_constant) * 0.95 * initial_voltage,
    (time, 3 * rc_time_constant, 3 * rc_time_constant + 0.001),
    label = 'time = 3 * Tau',
    line_color='green',
    show=False,
    )
UC.append(tau3_line[0])

voltageFull = plot(
    initial_voltage,
    (time, 0, 8 * rc_time_constant),
    line_color='red',      
    label='U0',
    backend=MatplotlibBackend,
    show=False)
UC.append(voltageFull[0])

UC.show()
'''