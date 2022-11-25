#!/usr/bin/env python3

# Description
## You throw the grenade to enemies. Which angle of throwing should you choose to hit the farest enemy?

from math import pi
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import symbols, Eq, solve, Function
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as movement_law
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector

flight_time = symbols('flight_time')
throwing_velocity = symbols('throwing_velocity')
throwing_angle = symbols('throwing_angle')

time_argument = symbols('time_argument')
gravitational_acceleration = symbols('gravitational_acceleration')

## Choose coordinates: object starts its flight from (0, 0), vertical axis is upwards, horisontal axis is toward object's destination.
## So the angle between X-axis and initial speed is throwing_angle, and the angle between Y-axis and gravitational acceleration is pi / 2 - throwing_angle.
## Both X and Y projections of movement are constant acceleration movements. 
## The vertical initial speed is projection of throwing_velocity to Y, acceleration is projection of gravitational acceleration to Y and it is -g.
## The horizontal initial speed is projection of throwing_velocity to X, acceleration is 0.
## Flight ends when y == 0 again.

initial_vertical_velocity = solve(projector.law, projector.projection, dict=True)[0][projector.projection].subs(
    {
        projector.vector_length: throwing_velocity, projector.vector_angle: pi / 2 - throwing_angle
    })

initial_horizontal_velocity = solve(projector.law, projector.projection, dict=True)[0][projector.projection].subs(
    {
        projector.vector_length: throwing_velocity, projector.vector_angle: throwing_angle
    })

x_function = symbols('x_function', cls=Function)
x_function = movement_law.law.rhs.subs(
    {        
        movement_law.initial_velocity: initial_horizontal_velocity, 
        movement_law.movement_time: time_argument, 
        movement_law.constant_acceleration: 0
    })

y_function = symbols('y_function', cls=Function)
y_function = movement_law.law.rhs.subs(
    {         
        movement_law.movement_time: time_argument, 
        movement_law.constant_acceleration: -gravitational_acceleration        
    })

end_of_flight = Eq(0, y_function)

# there are 2 points of the trajectory with h = 0, the first one is the start point, and the second one is at destination
# We need the second one
flight_time = solve(end_of_flight, time_argument, dict=True)[1][time_argument]

flight_distance = x_function.subs({time_argument: flight_time})
print(flight_distance)

# Let's find needed angle analytically.

'''
Len = plot(
    flight_distance,
    (throwing_angle, 0, pi / 2),
    line_color='blue',    
    title='Length',    
    label = 'Length',    
    legend=True,    
    annotations = {},
    backend=MatplotlibBackend,    
    show=False)

tau_line = plot(
    1000 * (time_argument - pi / 4),
    (time_argument, pi / 4, pi / 4 + 0.001),
    label = 'angle = pi / 4',
    line_color='yellow',
    show=False,
)
Len.append(tau_line[0])
Len.show()
'''