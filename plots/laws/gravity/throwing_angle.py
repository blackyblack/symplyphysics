#!/usr/bin/env python3

# Description
## You throw the grenade to enemies. Which angle of throwing should you choose to hit the farest enemy?

from sympy import diff, symbols, Eq, solve, simplify, pi
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as movement_law
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector

flight_time = symbols("flight_time")
throwing_velocity = symbols("throwing_velocity")
# hint sympy that we are not interested in negative throwing angles
throwing_angle = symbols("throwing_angle", positive=True)
time_argument = symbols("time_argument")
gravitational_acceleration = symbols("gravitational_acceleration")

## Choose coordinates: object starts its flight from (0, 0), vertical axis is upwards, horisontal axis is toward object's destination.
## So the angle between X-axis and initial speed is throwing_angle, and the angle between Y-axis and gravitational acceleration is pi / 2 - throwing_angle.
## Both X and Y projections of movement are constant acceleration movements. 
## The vertical initial speed is projection of throwing_velocity to Y, acceleration is projection of gravitational acceleration to Y and it is -g
## (-g is the projection of free fall acceleration which is downwards to the Y axis which is upwards, so angle between acceleration and Axis is 180°)
## The horizontal initial speed is projection of throwing_velocity to X, acceleration is 0.
## Flight ends when y == 0 again.

initial_horizontal_velocity = solve(projector.law, projector.projection, dict=True)[0][projector.projection].subs({
    projector.vector_length: throwing_velocity, projector.vector_angle: throwing_angle})
print(f"Initial horizontal velocity: {initial_horizontal_velocity}")

# the angle between initial velocity and Y axis is (pi/2 - throwing angle) (that's because of pi/2 angle between X and Y axis)
initial_vertical_velocity = solve(projector.law, projector.projection, dict=True)[0][projector.projection].subs({
    projector.vector_length: throwing_velocity, projector.vector_angle: pi / 2 - throwing_angle})
print(f"Initial vertical velocity: {initial_vertical_velocity}")

horizontal_movement = movement_law.law.rhs.subs({
    movement_law.initial_velocity: initial_horizontal_velocity, 
    movement_law.movement_time: time_argument, 
    movement_law.constant_acceleration: 0})
print(f"x(t): {horizontal_movement}")

vertical_movement = movement_law.law.rhs.subs({
    movement_law.initial_velocity: initial_vertical_velocity,
    movement_law.movement_time: time_argument, 
    movement_law.constant_acceleration: -gravitational_acceleration})
print(f"y(t): {vertical_movement}")

end_of_flight = Eq(0, vertical_movement)
# there are 2 points of the trajectory with h = 0, the first one is the start point, and the second one is at destination
# We need the second one
flight_time = solve(end_of_flight, time_argument, dict=True)[1][time_argument]
print(f"Flight time: {flight_time}")

flight_distance = horizontal_movement.subs({time_argument: flight_time})
print(f"Flight distance: {flight_distance}")

# for plotting purposes we don't care about scale and may substitute constants with 1
flight_distance_plotted = flight_distance.subs({
    throwing_velocity: 1,
    gravitational_acceleration: 1})
flight_distance_plotted = simplify(flight_distance_plotted)
print(f"Plotted flight distance: {flight_distance_plotted}")

# Let's find needed angle analytically.
# We are to find maximum of flight_distance(throwing_angle) function with angle in (0, pi/2) interval.
distance_diff = diff(flight_distance, throwing_angle)
max_law = Eq(0, distance_diff)
max_angle = solve(max_law, throwing_angle, dict=True)[0][throwing_angle]
print(f"Angle to achieve maximum distance: {max_angle}")

p0 = plot(
    flight_distance_plotted,
    (throwing_angle, 0, pi / 2),
    line_color="blue",
    title="Throwing distance depending on angle",
    label = "distance",
    xlabel="throwing angle, radians",
    ylabel="distance",
    legend=True,
    annotations = {},
    backend=MatplotlibBackend,
    show=False)

peak_line = plot(
    1000 * (throwing_angle - max_angle),
    (throwing_angle, max_angle, max_angle + 0.001),
    label = "angle = 45°",
    line_color="green",
    show=False,
)
p0.append(peak_line[0])
p0.show()
