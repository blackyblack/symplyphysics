#!/usr/bin/env python3
from symplyphysics.acceleration_from_speed_difference import module as acceleration_definition

v0 = 1
v1 = 20
t = 5
print("Formula is:\n{}".format(acceleration_definition.print()))
print("Acceleration = {}; for initial velocity = {}, terminal velocity = {}, time period = {}"
    .format(acceleration_definition.calculate_linear_acceleration(v0, v1, t), v0, v1, t))