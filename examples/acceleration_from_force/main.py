#!/usr/bin/env python3
from symplyphysics.laws.dynamics import acceleration_from_force as newton_law2

m = 1
a = 3
print("Formula is:\n{}".format(newton_law2.print()))
print("Force = {}; for mass = {}, acceleration = {}".format(newton_law2.calculate_force(m, a), m, a))