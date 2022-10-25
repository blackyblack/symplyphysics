#!/usr/bin/env python3

# Description
## You throw the grenade to enemies. Which angle of throwing should you choose to hit the farest enemy?
## Flying velocity V0 vector is the sum of two orthogonal - vertical and horizontal speed vectors V_vert and V_hor.
## If alfa is throwing angle, V_vert = V0 * sin(alpha) and V_hor = V0 * cos(alpha)
## Vertical speed is accelerated by the gravitational acceleration.
## Horizontal speed is constant.

from math import pi
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

print(flight_time) #should be 2 solvations, we need one which is not zero
##kostil
flight_time = 2 * V_vertical / g

Length = symbols('Length', cls = Function)
Length = (V_horizontal * flight_time).subs({V0: 1})

print(Length)

Len = plot(
    Length,
    (alpha, 0, pi / 2),
    line_color='blue',    
    title='Length',    
    label = 'Length',    
    legend=True,    
    annotations = {},
    backend=MatplotlibBackend,    
    show=False)

tau_line = plot(
    100 * (time - pi / 4),
    (time, pi / 4, pi / 4 + 0.001),
    label = 'angle = pi / 4',
    line_color='yellow',
    show=False,
)
Len.append(tau_line[0])
Len.show()
