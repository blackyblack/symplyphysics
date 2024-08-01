#!/usr/bin/env python3
from sympy import sin
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.definitions import acceleration_is_speed_derivative as acceleration

print(f"Formula is:\n{print_expression(acceleration.definition)}")

VelocityFunction = sin
applied_law = acceleration.definition.subs(acceleration.speed, VelocityFunction)
dsolved = applied_law.doit()

print(f"Velocity function is:\n{print_expression(VelocityFunction(acceleration.time))}")
print(f"Acceleration function is:\n{print_expression(dsolved)}")

p1 = plot(VelocityFunction(acceleration.time), (acceleration.time, 0, 10),
    line_color="blue",
    title="Acceleration(time), Velocity(time)",
    label="Velocity(time)",
    xlabel="time",
    ylabel="f(time)",
    legend=True,
    backend=MatplotlibBackend,
    show=False)

p2 = plot(dsolved.rhs, (acceleration.time, 0, 10),
    line_color="red",
    label="Acceleration(time)",
    backend=MatplotlibBackend,
    show=False)

p1.append(p2[0])
p1.show()
