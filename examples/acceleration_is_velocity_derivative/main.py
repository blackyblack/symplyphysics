#!/usr/bin/env python3
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration

v0 = 1
v1 = 20
t = 5
print("Formula is:\n{}".format(acceleration.print()))
print("Acceleration = {}; for initial velocity = {}, terminal velocity = {}, time period = {}"
    .format(acceleration.calculate_linear_acceleration(v0, v1, t), v0, v1, t))