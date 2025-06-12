#!/usr/bin/env python3

from collections import namedtuple
from sympy import Symbol, plot, sqrt, Expr
from symplyphysics.quantities import speed_of_light
from symplyphysics.core.convert import evaluate_expression
from symplyphysics.laws.relativistic.vector import force_acceleration_relation as force_law
from symplyphysics.definitions import lorentz_factor as lorentz_factor_def

from symplyphysics.core.experimental.vectors import VectorNorm, VectorDot
from symplyphysics.core.experimental.coordinate_systems import CoordinateVector, CARTESIAN

Datum = namedtuple("Datum", "label acceleration")

data_ = [
    Datum(
    label=r"$a_{\perp} = 0$",
    acceleration=CoordinateVector([1, 0, 0], CARTESIAN),
    ),
    Datum(
    label=r"$a_{\perp} \ne 0, a_{\parallel} \ne 0$",
    acceleration=CoordinateVector([1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)], CARTESIAN),
    ),
    Datum(
    label=r"$a_{\parallel} = 0$",
    acceleration=CoordinateVector([0, 1 / sqrt(2), 1 / sqrt(2)], CARTESIAN),
    ),
]

reduced_speed = Symbol("reduced_speed", nonnegative=True)
velocity = CoordinateVector([reduced_speed * speed_of_light, 0, 0], CARTESIAN)

base_plot = plot(
    title="Force as a function of velocity for different acceleration configurations",
    xlabel=r"$\beta = \frac{v}{c}$",
    ylabel=r"$\frac{F}{m_0}, \left(\frac{\text{m}}{\text{s}}\right)^2$",
    legend=True,
    show=False,
)

v_dot_v = VectorDot(velocity, velocity)


def split_acceleration(acceleration: Expr) -> tuple[Expr, Expr]:
    tangential = VectorDot(acceleration, velocity) / v_dot_v * velocity

    normal = acceleration - tangential

    return tangential, normal


lorentz_factor = lorentz_factor_def.definition.rhs.subs(
    lorentz_factor_def.speed,
    reduced_speed * speed_of_light,
)

for datum_ in data_:
    tangential_acceleration, normal_acceleration = split_acceleration(datum_.acceleration)

    vector = force_law.law.rhs.subs({
        force_law.rest_mass: 1,
        force_law.tangential_acceleration: tangential_acceleration,
        force_law.normal_acceleration: normal_acceleration,
        force_law.lorentz_factor: lorentz_factor,
    })
    vector = CoordinateVector.from_expr(vector)

    expr = VectorNorm(vector)
    expr_value = evaluate_expression(expr)
    sub_plot = plot(
        expr_value,
        (reduced_speed, 0, 0.97),
        label=datum_.label,
        show=False,
    )
    base_plot.extend(sub_plot)

base_plot.show()
