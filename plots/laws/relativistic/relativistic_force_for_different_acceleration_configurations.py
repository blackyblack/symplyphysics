#!/usr/bin/env python3

from collections import namedtuple
from sympy import Symbol, plot, sqrt
from sympy.physics.units import speed_of_light
from symplyphysics import Vector, vector_magnitude
from symplyphysics.laws.relativistic.vector import force_acceleration_relation

Datum = namedtuple("Datum", "label acceleration")

data_ = [
    Datum(
        label=r"$a_{\perp} = 0$",
        acceleration=Vector([1, 0, 0]),
    ),
    Datum(
        label=r"$a_{\perp} \ne 0, a_{\parallel} \ne 0$",
        acceleration=Vector([1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)]),
    ),
    Datum(
        label=r"$a_{\parallel} = 0$",
        acceleration=Vector([0, 1 / sqrt(2), 1 / sqrt(2)]),
    ),
]

reduced_velocity = Symbol("reduced_velocity", nonnegative=True)
velocity = Vector([reduced_velocity * speed_of_light, 0, 0])

base_plot = plot(
    title="Force as a function of velocity for different acceleration configurations",
    xlabel=r"$\beta = \frac{v}{c}$",
    ylabel=r"$\frac{F}{m_0}, \left(\frac{\text{m}}{\text{s}}\right)^2$",
    legend=True,
    show=False,
)

for datum_ in data_:
    vector = force_acceleration_relation.force_law(datum_.acceleration, velocity)
    expr = vector_magnitude(vector).subs(force_acceleration_relation.rest_mass, 1)

    sub_plot = plot(
        expr,
        (reduced_velocity, 0, 0.97),
        label=datum_.label,
        show=False,
    )
    base_plot.extend(sub_plot)

base_plot.show()
