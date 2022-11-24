#!/usr/bin/env python3

# Description
## You throw the grenade to enemies. Which angle of throwing should you choose to hit the farest enemy?

from math import pi
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import symbols, Eq, pretty, solve, Function
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as movement_law
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector

flight_time = symbols('flight_time')
throwing_velocity = symbols('throwing_velocity')
alpha = symbols('alpha')
Distance = symbols('Distance')

velocity_vertical = projector.calculate_projection(throwing_velocity, pi / 2 - alpha)
velocity_horizontal = projector.calculate_projection(throwing_velocity, alpha)
time = symbols('time')
gravitational_constant = symbols('gravitational_constant')

## Height function of time is h = V_vertical * t - g * t**2 /2
## Flight ends when h == 0

end_of_flight = Eq(0, velocity_vertical * time - gravitational_constant * time**2 /2)
flight_time = solve(end_of_flight, time, dict=True)[1][time] 
# there are 2 points of the trajectory with h = 0, the first one is the start point, and the second one is at destination
# We need the second one

Length = (velocity_horizontal * flight_time).subs({throwing_velocity: 1, gravitational_constant: 1})

# Let's find needed angle analytically.


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
    1000 * (time - pi / 4),
    (time, pi / 4, pi / 4 + 0.001),
    label = 'angle = pi / 4',
    line_color='yellow',
    show=False,
)
Len.append(tau_line[0])
Len.show()
