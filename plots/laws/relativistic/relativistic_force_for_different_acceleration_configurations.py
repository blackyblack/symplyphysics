#!/usr/bin/env python3

from collections import namedtuple
from sympy import Symbol, plot
from sympy.physics.units import speed_of_light
from symplyphysics import Vector, vector_magnitude
from symplyphysics.laws.relativistic.vector import force_via_acceleration

Datum = namedtuple("Datum", "label acceleration")

data_ = [
    Datum(
        label=r"$a_{\perp} = 0$",
        acceleration=Vector([10, 0, 0]),
    ),
    Datum(
        label=r"$a_{\parallel} = 0$",
        acceleration=Vector([0, 10, -20]),
    ),
    Datum(
        label=r"$a_{\perp} \ne 0, a_{\parallel} \ne 0$",
        acceleration=Vector([10, 10, -20]),
    ),
]

reduced_velocity = Symbol("reduced_velocity", nonnegative=True)
velocity = Vector([reduced_velocity * speed_of_light, 0, 0])

base_plot = plot(
    title="Force as a function of velocity for different acceleration configurations",
    xlabel=r"$\beta = \frac{v}{c}$",
    ylabel=r"$\frac{F}{m_0}, \frac{\text{m}^2}{\text{s}^2}$",
    legend=True,
    show=False,
)

for datum_ in data_:
    vector = force_via_acceleration.force_law(datum_.acceleration, velocity)
    expr = vector_magnitude(vector).subs(force_via_acceleration.rest_mass, 1)

    sub_plot = plot(
        expr,
        (reduced_velocity, 0, 0.96),
        label=datum_.label,
        show=False,
    )
    base_plot.extend(sub_plot)

base_plot.show()
