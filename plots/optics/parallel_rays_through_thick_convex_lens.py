#!/usr/bin/env python3
"""
Plot the passage of rays parallel to the optical axis through a thick convex lens.
"""

from __future__ import annotations
from typing import Any
from math import isnan, isinf
from argparse import ArgumentParser
import sys
from sympy import (ImmutableMatrix as Matrix, symbols as sym_symbols, Eq, solve, Expr, sqrt,
    rot_axis3, sin, pi, asin, Point2D, sign)
from sympy.plotting import plot, plot_parametric
from symplyphysics.laws.optics import refraction_angle_from_environments as snells_law
from symplyphysics.core.geometry.line import two_point_function

# PART 1. Parsing command line arguments


def make_parser() -> ArgumentParser:
    parser = ArgumentParser()
    parser.add_argument(
        "-r1",
        action="store",
        default=2.0,
        dest="r1",
        help="radius of curvature of the left side of the lens",
    )
    parser.add_argument(
        "-r2",
        action="store",
        default=2.0,
        dest="r2",
        help="radius of curvature of the right side of the lens",
    )
    parser.add_argument(
        "-n1",
        action="store",
        default=1.0,
        dest="n1",
        help="refractive index of the surrounding media, default: air (n=1)",
    )
    parser.add_argument(
        "-n2",
        action="store",
        default=1.5,
        dest="n2",
        help="refraction index of the lens material, default: crown glass (n=1.5)",
    )
    return parser


def positive_finite_float(s: str) -> float:
    v = float(s)
    if isnan(v):
        raise ValueError(f"Expect a non-NaN input, got '{s}'")
    if isinf(v):
        raise ValueError(f"Expect a finite input, got '{s}'")
    if v <= 0:
        raise ValueError(f"Expect a positive input, got '{s}'")
    return v


args = make_parser().parse_args(sys.argv[1:])

# Radius of curvature of the **left** side of the lens
r1 = positive_finite_float(args.r1)

# Radius of curvature of the **right** side of the lens
r2 = positive_finite_float(args.r2)

# Refractive index of the medium surrounding the lens
n1 = positive_finite_float(args.n1)

# Refractive index of the material of the lens
n2 = positive_finite_float(args.n2)

if n1 >= n2:
    MESSAGE = "Expect the refraction index of the surrounding medium to be smaller than that of the lens."
    raise ValueError(MESSAGE)

# PART 2. Calculations


def make_point(x_: Any, y_: Any) -> Matrix:
    """Embeds a 2-dimensional point `(x, y)` into the 3-dimensional plane `z = 0`."""

    return Matrix([x_, y_, 0])


# Radius of the lens
r = 1.0

# Thickness of the lens
d = (r1 - sqrt(r1**2 - r**2)) + (r2 - sqrt(r2**2 - r**2))

# Effective focus length of the lens
# TODO Use lensmaker equation to calculate effective focus length:
f = 1 / ((n2 - n1) * (1 / r1 + 1 / r2 - (n2 - n1) * d / (n2 * r1 * r2)))

# Focus length of an ideal thin lens with the same radii of curvature
f_thin = 1 / ((n2 - n1) * (1 / r1 + 1 / r2))

# Position vector of the center of the circle making up the **left** side of the lens
c1 = make_point(+sqrt(r1**2 - r**2), 0)

# Position vector of the center of the circle making up the **right** side of the lens
c2 = make_point(-sqrt(r2**2 - r**2), 0)

x, y = sym_symbols("x y", real=True)
p = make_point(x, y)


def sq_norm(v: Matrix) -> Expr:
    """Returns the square of the norm of `v`."""
    return v.dot(v)


# Equation of the **left** side of the lens
eqn1 = Eq(sq_norm(p - c1), r1**2)

# Equation of the **right** side of the lens
eqn2 = Eq(sq_norm(p - c2), r2**2)


def normalize(v: Matrix) -> Expr:
    """Returns the unit direction vector co-directed with the given vector."""
    return v / sqrt(sq_norm(v))


def sin_between_unit_vectors(u1: Matrix, u2: Matrix) -> Expr:
    """Assumes unit input vectors with a zero third component."""
    return abs(u1.cross(u2)[2])


def rotate(theta: Expr, v: Matrix) -> Matrix:
    """
    Rotates vector `v` clockwise at angle `theta` if `theta > 0`, and anti-clockwise otherwise.

    >>> assert rotate(pi / 2, make_point(1, 0)) == make_point(0, 1)
    >>> assert rotate(-pi / 2, make_point(1, 1)) == make_point(1, -1)
    """

    return rot_axis3(-theta) * v


def calculate(y_in: float) -> tuple[Matrix, Matrix, Matrix] | None:
    """Calculates intersections with the lens and with the optical axis `y = 0`"""

    # Point of intersection of the incoming ray with the left side of the lens
    (x_in,) = (x for x in solve(eqn1.subs(y, y_in), x) if x < 0)
    p_in = make_point(x_in, y_in)

    # Unit normal on the left surface of the lens at the intersection `p_in`
    n_in = normalize(p_in - c1)

    # Unit vector of the incoming ray
    u_in1 = make_point(1, 0)
    sin_in1 = sin_between_unit_vectors(n_in, u_in1)

    snells_eqn_in = snells_law.law.subs({
        snells_law.incidence_refractive_index: n1,
        sin(snells_law.incidence_angle): sin_in1,
        snells_law.resulting_refractive_index: n2,
    })
    (theta_in2,) = (
        theta for theta in solve(snells_eqn_in, snells_law.refraction_angle) if theta < pi / 2)

    # Unit vector of the ray traveling inside the lens
    u2 = rotate(sign(y_in) * theta_in2, -n_in)
    ray2_eqn = Eq((p - p_in).cross(u2)[2], 0)

    # Point of intersection of the ray in the lens and the right side of the lens.
    (p_out,) = (make_point(sol[x], sol[y])
        for sol in solve((eqn2, ray2_eqn), (x, y), dict=True)
        if sol[x] > 0)

    # Unit normal on the right surface of the lens at the intersection `p_out`
    n_out = normalize(p_out - c2)
    sin_out2 = sin_between_unit_vectors(n_out, u2)

    snells_eqn_out = snells_law.law.subs({
        snells_law.incidence_refractive_index: n2,
        sin(snells_law.incidence_angle): sin_out2,
        snells_law.resulting_refractive_index: n1,
    })
    solved = [
        theta for theta in solve(snells_eqn_out, snells_law.refraction_angle)
        if abs(theta) <= pi / 2
    ]
    if not solved:
        # Case of total internal reflection
        return None
    (theta_out1,) = solved

    # Unit vector of the outgoing ray
    u_out1 = rotate(-sign(y_in) * theta_out1, n_out)

    ray_out1_eqn = Eq((p - p_out).cross(u_out1)[2], 0)

    # Point of intersection of the outgoing ray with the optical axis
    x_f = solve(ray_out1_eqn.subs(y, 0), x)[0]
    p_f = make_point(x_f, 0)

    return p_in, p_out, p_f


# PART 3. Plotting


def to_point2d(m: Matrix) -> Point2D:
    return Point2D(m[0], m[1])


def make_line_eqn(v1: Matrix, v2: Matrix) -> Expr:
    return two_point_function(to_point2d(v1), to_point2d(v2), x)


def display() -> None:
    y_ins = (-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9)

    data: list[tuple[Matrix, Matrix, Matrix]] = []
    for y_in in y_ins:
        t = calculate(y_in)
        if not t:
            continue
        data.append(t)

    maxval = max(max((t[2][0] for t in data), default=1), 1) * 1.1

    base_plot = plot(
        title="Parallel rays passing through a thick convex lens",
        xlabel="",
        ylabel="",
        ylim=(-maxval, maxval),
        xlim=(-maxval, maxval),
        size=(8, 8),
        legend=True,
        show=False,
    )

    subplot = plot_parametric(
        (f, y),
        (y, -maxval, maxval),
        line_color="pink",
        label="effective focus plane",
        show=False,
    )
    base_plot.extend(subplot)

    phi = sym_symbols("phi", real=True)
    phi1 = pi - asin(r / r1)
    phi2: Expr = asin(r / r2)

    eqn1_parametric = rotate(phi, make_point(r1, 0)) + c1
    subplot = plot_parametric(
        eqn1_parametric[:-1],
        (phi, phi1, 2 * pi - phi1),
        line_color="blue",
        label="",
        show=False,
    )
    base_plot.extend(subplot)

    eqn2_parametric = rotate(phi, make_point(r2, 0)) + c2
    subplot = plot_parametric(
        eqn2_parametric[:-1],
        (phi, -phi2, phi2),
        line_color="blue",
        label="",
        show=False,
    )
    base_plot.extend(subplot)

    for p_in, p_out, p_f, *_ in data:
        p_left = make_point(-maxval, p_in[1])
        subplot = plot(
            make_line_eqn(p_left, p_in),
            (x, -maxval, p_in[0]),
            line_color="green",
            label="",
            show=False,
        )
        base_plot.extend(subplot)

        subplot = plot(
            make_line_eqn(p_in, p_out),
            (x, p_in[0], p_out[0]),
            line_color="green",
            label="",
            show=False,
        )
        base_plot.extend(subplot)

        subplot = plot(
            make_line_eqn(p_out, p_f),
            (x, p_out[0], p_f[0]),
            line_color="green",
            label="",
            show=False,
        )
        base_plot.extend(subplot)

    base_plot.show()


display()
