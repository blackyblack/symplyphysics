#!/usr/bin/env python3
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration
from sympy import sin, solve, dsolve, pretty
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend

print("Formula is:\n{}".format(acceleration.print()))

velocity_function = sin
applied_law = acceleration.law.subs(acceleration.velocity_function, velocity_function)
dsolved = dsolve(applied_law, acceleration.acceleration(acceleration.time))

print("Velocity function is:\n{}".format(pretty(velocity_function(acceleration.time), use_unicode=False)))
print("Acceleration function is:\n{}".format(pretty(dsolved, use_unicode=False)))

p1 = plot(
    velocity_function(acceleration.time),
    (acceleration.time, 0, 10),
    line_color='blue',
    title='Acceleration(time), Velocity(time)',
    label='Velocity(time)',
    legend=True,
    backend=MatplotlibBackend,
    show=False)

p2 = plot(
    dsolved.rhs,
    (acceleration.time, 0, 10),
    line_color='red',
    label='Acceleration(time)',
    backend=MatplotlibBackend,
    show=False)

p1.append(p2[0])
p1.show()
