#!/usr/bin/env python3

from collections import namedtuple
from sympy import solve, symbols, Eq
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from symplyphysics.quantities import speed_of_light
from symplyphysics.laws.relativistic import coordinate_conversion_velocity_constant as transform_law

reduced_speed = symbols("beta")
speed = transform_law.proper_frame_speed_in_lab_frame

reduced_speed_eqn = Eq(reduced_speed, speed / speed_of_light)

Datum = namedtuple("Datum", "x t label")

data_ = [
    Datum(x=0.5, t=-1 / speed_of_light, label=r"$x = \frac{1}{2}, ct = -1$"),
    Datum(x=-1, t=-1 / speed_of_light, label=r"$x = -1, ct = -1$"),
    Datum(x=1, t=2 / speed_of_light, label=r"$x = 1, ct = 2$"),
    Datum(x=-2, t=-1 / speed_of_light, label=r"$x = -2, ct = -1$"),
]

proper_coordinate_expr = solve(
    (transform_law.law, reduced_speed_eqn),
    (transform_law.position_in_proper_frame, speed),
    dict=True,
)[0][transform_law.position_in_proper_frame]

print(print_expression(proper_coordinate_expr))

base_plot = plot(
    title=r"Proper coordinate $x'$ as a function of frame speed",
    xlabel=r"$\beta = \frac{v}{c}$",
    ylabel=r"$x'$",
    backend=MatplotlibBackend,
    legend=True,
    annotation=False,
    show=False,
)

for datum_ in data_:
    expr = proper_coordinate_expr.subs({
        transform_law.position_in_lab_frame: datum_.x,
        transform_law.time_in_lab_frame: datum_.t,
    })
    sub_plot = plot(
        expr,
        (reduced_speed, 0, 0.98),
        label=datum_.label,
        show=False,
    )
    base_plot.extend(sub_plot)

base_plot.show()
