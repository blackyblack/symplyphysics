#!/usr/bin/env python3
from symplyphysics.acceleration_from_speed_difference import module as acceleration_definition
from sympy import sin, solve, dsolve, pretty
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend

velocity_function = sin

print("Formula is:\n{}".format(acceleration_definition.print()))

applied_law = acceleration_definition.law.subs(acceleration_definition.velocity_function_, velocity_function)
dsolved = dsolve(applied_law, acceleration_definition.acceleration_(acceleration_definition.time_))

print("Acceleration function is:\n{}".format(pretty(dsolved, use_unicode=False)))

p1 = plot(
    velocity_function(acceleration_definition.time_),
    (acceleration_definition.time_, 0, 10),
    line_color='blue',
    title='Acceleration(time), Velocity(time)',
    label='Velocity(time)',
    legend=True,
    backend=MatplotlibBackend,
    show=False)

p2 = plot(
    dsolved.rhs,
    (acceleration_definition.time_, 0, 10),
    line_color='red',
    label='Acceleration(time)',
    backend=MatplotlibBackend,
    show=False)

p1.append(p2[0])
p1.show()