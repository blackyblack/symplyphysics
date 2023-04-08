#!/usr/bin/env python3
from sympy import sin, pretty
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics.core.symbols.symbols import to_printable
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration

print("Formula is:\n{}".format(acceleration.print()))

velocity_function = sin
applied_law = acceleration.definition.subs(acceleration.velocity_function, velocity_function)
dsolved = applied_law.doit()

print("Velocity function is:\n{}".format(pretty(to_printable(velocity_function(acceleration.time), acceleration.symbols), use_unicode=False)))
print("Acceleration function is:\n{}".format(pretty(to_printable(dsolved, acceleration.symbols), use_unicode=False)))

p1 = plot(
    velocity_function(acceleration.time),
    (acceleration.time, 0, 10),
    line_color='blue',
    title='Acceleration(time), Velocity(time)',
    label='Velocity(time)',
    xlabel='time',
    ylabel='f(time)',
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
