#!/usr/bin/env python3
from sympy import sin
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration

print("Formula is:\n{}".format(acceleration.print()))

velocity_function = sin
applied_law = acceleration.definition.subs(acceleration.velocity, velocity_function)
dsolved = applied_law.doit()

print("Velocity function is:\n{}".format(print_expression(velocity_function(acceleration.time))))
print("Acceleration function is:\n{}".format(print_expression(dsolved)))

p1 = plot(
    velocity_function(acceleration.time),
    (acceleration.time, 0, 10),
    line_color="blue",
    title="Acceleration(time), Velocity(time)",
    label="Velocity(time)",
    xlabel="time",
    ylabel="f(time)",
    legend=True,
    backend=MatplotlibBackend,
    show=False)

p2 = plot(
    dsolved.rhs,
    (acceleration.time, 0, 10),
    line_color="red",
    label="Acceleration(time)",
    backend=MatplotlibBackend,
    show=False)

p1.append(p2[0])
p1.show()
