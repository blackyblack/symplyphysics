#!/usr/bin/env python3
"""
Plot the passage of parallel rays through a thick convex lens.
"""

from __future__ import annotations
from typing import Any
from math import isnan
from argparse import ArgumentParser
import sys
from sympy import (ImmutableMatrix as Matrix, symbols as sym_symbols, Eq, solve, Expr, sqrt,
    rot_axis3, sin, pi, asin, Point2D, sign)
from sympy.plotting import plot, plot_parametric
from symplyphysics.laws.optics import refraction_angle_from_environments as snells_law
from symplyphysics.core.geometry.line import two_point_function

# PART 1. Parsing command line arguments


def positive_float(s: str) -> float:
    f = float(s)
    if isnan(f):
        raise ValueError(f"Expect a non-NaN input, got '{s}'")
    if f <= 0:
        raise ValueError(f"Expect positive input, got '{f}'")
    return f


def make_parser() -> ArgumentParser:
    parser = ArgumentParser()
    parser.add_argument(
        "-r1",
        action="store",
        default=10.0,
        type=positive_float,
        dest="r1",
        help="radius of curvature of the left side of the lens",
    )
    parser.add_argument(
        "-r2",
        action="store",
        default=10.0,
        type=positive_float,
        dest="r2",
        help="radius of curvature of the right side of the lens",
    )
    parser.add_argument(
        "-n1",
        action="store",
        default=1.0,
        type=positive_float,
        dest="n1",
        help="refractive index of the surrounding media (default: air)",
    )
    parser.add_argument(
        "-n2",
        action="store",
        default=1.5,
        type=positive_float,
        dest="n2",
        help="refraction index of the lens material (default: crown glass)",
    )
    return parser


args = make_parser().parse_args(sys.argv[1:])
r1, r2, n1, n2 = args.r1, args.r2, args.n1, args.n2

# PART 2. Calculations


def make_point(x_: Any, y_: Any) -> Matrix:
    return Matrix([x_, y_, 0])


r = 1.0
# TODO: add calculations of lens thickness `d` and effective focus length `f`
# d = (r1 - sqrt(r1**2 - r**2)) + (r2 - sqrt(r2**2 - r**2))
# 1 / f = (n2 - n1) * (1 / r1 + 1 / r2 - (n2 - n1) * d / (n2 * r1 * r2))

c1 = make_point(+sqrt(r1**2 - r**2), 0)
c2 = make_point(-sqrt(r2**2 - r**2), 0)
x, y = sym_symbols("x y", real=True)
p = make_point(x, y)


def sq_norm(m: Matrix) -> Expr:
    return m.dot(m)


def normalize(m: Matrix) -> Expr:
    return m / sqrt(sq_norm(m))


def sin_between_unit_vectors(u1: Matrix, u2: Matrix) -> Expr:
    """Assumes unit input vectors with a zero third component."""
    return abs(u1.cross(u2)[2])


def rotate(theta: Expr, v: Matrix) -> Matrix:
    return rot_axis3(theta) * v


eqn1 = Eq(sq_norm(p - c1), r1**2)
eqn2 = Eq(sq_norm(p - c2), r2**2)


def calculate(y_in: float) -> tuple[Matrix, Matrix, Matrix]:
    # TODO: add variable descriptions
    # NOTE: check calculations

    (x_in,) = (x for x in solve(eqn1.subs(y, y_in), x) if x < 0)
    p_in = make_point(x_in, y_in)

    n_in = normalize(p_in - c1)
    u_in1 = make_point(1, 0)
    sin_in1 = sin_between_unit_vectors(n_in, u_in1)

    snells_eqn_in = snells_law.law.subs({
        snells_law.incidence_refractive_index: n1,
        sin(snells_law.incidence_angle): sin_in1,
        snells_law.resulting_refractive_index: n2,
    })
    (theta_in2,) = (
        theta for theta in solve(snells_eqn_in, snells_law.refraction_angle) if theta < pi / 2)

    u2 = rotate(pi + theta_in2, n_in)
    ray2_eqn = Eq((p - p_in).cross(u2)[2], 0)

    (p_out,) = (p for p in solve((eqn2, ray2_eqn), (x, y), dict=True) if p[x] > 0)
    p_out = make_point(p_out[x], p_out[y])

    n_out = normalize(p_out - c2)
    sin_out2 = sin_between_unit_vectors(n_out, u2)

    snells_eqn_out = snells_law.law.subs({
        snells_law.incidence_refractive_index: n2,
        sin(snells_law.incidence_angle): sin_out2,
        snells_law.resulting_refractive_index: n1,
    })
    (theta_out1,) = (theta for theta in solve(snells_eqn_out, snells_law.refraction_angle)
        if abs(theta) < pi / 2)

    u_out1 = rotate(sign(y_in) * theta_out1, n_out)

    ray_out1_eqn = Eq((p - p_out).cross(u_out1)[2], 0)
    x_f = solve(ray_out1_eqn.subs(y, 0), x)[0]
    p_f = make_point(x_f, 0)

    return p_in, p_out, p_f, u_out1


# PART 3. Plotting


def to_point2d(m: Matrix) -> Point2D:
    return Point2D(m[0], m[1])


def make_line_eqn(v1: Matrix, v2: Matrix) -> Expr:
    return two_point_function(to_point2d(v1), to_point2d(v2), x)


def display() -> None:
    # FIXME: fix positive range of `y_in`
    y_ins = (-0.9, -0.6, -0.3, 0.3, 0.6, 0.9)
    data = tuple(map(calculate, y_ins))
    maxval = max(max(t[2][0] for t in data), 1.1)
    print(maxval)

    base_plot = plot(
        title="Parallel rays passing a thick convex lens",
        xlabel="",
        ylabel="",
        ylim=(-maxval, maxval),
        xlim=(-maxval, maxval),
        show=False,
    )

    phi = sym_symbols("phi", real=True)
    phi1 = pi - asin(r / r1)
    phi2: Expr = asin(r / r2)

    eqn1_parametric = rotate(phi, make_point(r1, 0)) + c1
    base_plot.extend(
        plot_parametric(
        eqn1_parametric[:-1],
        (phi, phi1, 2 * pi - phi1),
        show=False,
        line_color="blue",
        ))

    eqn2_parametric = rotate(phi, make_point(r2, 0)) + c2
    base_plot.extend(
        plot_parametric(
        eqn2_parametric[:-1],
        (phi, -phi2, phi2),
        show=False,
        line_color="blue",
        ))

    for p_in, p_out, p_f, *_ in data:
        p_left = make_point(-maxval, p_in[1])
        base_plot.extend(
            plot(
            make_line_eqn(p_left, p_in),
            (x, -maxval, p_in[0]),
            line_color="green",
            show=False,
            ))

        base_plot.extend(
            plot(
            make_line_eqn(p_in, p_out),
            (x, p_in[0], p_out[0]),
            line_color="green",
            show=False,
            ))

        base_plot.extend(
            plot(
            make_line_eqn(p_out, p_f),
            (x, p_out[0], p_f[0]),
            line_color="green",
            show=False,
            ))

    base_plot.show()


display()
