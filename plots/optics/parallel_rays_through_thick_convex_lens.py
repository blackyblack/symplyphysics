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
    rot_axis3, sin, pi, asin, Point2D, sign, im)
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
        dest="left_radius",
        help="radius of curvature of the left side of the lens",
    )
    parser.add_argument(
        "-r2",
        action="store",
        default=2.0,
        dest="right_radius",
        help="radius of curvature of the right side of the lens",
    )
    parser.add_argument(
        "-n1",
        action="store",
        default=1.0,
        dest="media_refractive_index",
        help="refractive index of the surrounding media, default: air (n=1)",
    )
    parser.add_argument(
        "-n2",
        action="store",
        default=1.5,
        dest="lens_refractive_index",
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
left_radius = positive_finite_float(args.left_radius)

# Radius of curvature of the **right** side of the lens
right_radius = positive_finite_float(args.right_radius)

# Refractive index of the medium surrounding the lens
media_refractive_index = positive_finite_float(args.media_refractive_index)

# Refractive index of the material of the lens
lens_refractive_index = positive_finite_float(args.lens_refractive_index)

if media_refractive_index >= lens_refractive_index:
    MESSAGE = "Expect the refraction index of the surrounding medium to be smaller than that of the lens."
    raise ValueError(MESSAGE)

# PART 2. Calculations


def make_point(x_: Any, y_: Any) -> Matrix:
    """Embeds a 2-dimensional point `(x, y)` into the 3-dimensional plane `z = 0`."""

    return Matrix([x_, y_, 0])


# Radius of the lens
lens_radius = 1.0

# Thickness of the lens
lens_thickness = (left_radius - sqrt(left_radius**2 - lens_radius**2)) + (right_radius -
    sqrt(right_radius**2 - lens_radius**2))

# Effective focus length of the lens
# TODO Use lensmaker equation to calculate effective focus length:
effective_focus_length = 1 / ((lens_refractive_index - media_refractive_index) *
    (1 / left_radius + 1 / right_radius -
    (lens_refractive_index - media_refractive_index) * lens_thickness /
    (lens_refractive_index * left_radius * right_radius)))

# Focus length of an ideal thin lens with the same radii of curvature
ideal_focus_length = 1 / ((lens_refractive_index - media_refractive_index) *
    (1 / left_radius + 1 / right_radius))

# Position vector of the center of the circle making up the **left** side of the lens
left_side_center = make_point(+sqrt(left_radius**2 - lens_radius**2), 0)

# Position vector of the center of the circle making up the **right** side of the lens
right_side_center = make_point(-sqrt(right_radius**2 - lens_radius**2), 0)

x, y = sym_symbols("x y", real=True)
p = make_point(x, y)


def sq_norm(v: Matrix) -> Expr:
    """Returns the square of the norm of `v`."""
    return v.dot(v)


# Equation of the **left** side of the lens
left_side_eqn = Eq(sq_norm(p - left_side_center), left_radius**2)

# Equation of the **right** side of the lens
right_side_eqn = Eq(sq_norm(p - right_side_center), right_radius**2)


def normalize(v: Matrix) -> Expr:
    """Returns the unit direction vector co-directed with the given vector."""
    return v / sqrt(sq_norm(v))


e_z = Matrix([0, 0, 1])


def sin_between_unit_vectors(u1: Matrix, u2: Matrix) -> Expr:
    """Assumes unit input vectors with a zero third component."""
    return abs(u1.cross(u2).dot(e_z))


def rotate(theta: Expr, v: Matrix) -> Matrix:
    """
    Rotates vector `v` anti-clockwise at angle `theta` if `theta > 0`, and clockwise otherwise.

    >>> assert rotate(pi / 2, make_point(1, 0)) == make_point(0, 1)
    >>> assert rotate(-pi / 2, make_point(1, 1)) == make_point(1, -1)
    """

    return rot_axis3(-theta) * v


(refraction_angle_expr,) = (
    sol for sol in solve(snells_law.law, snells_law.refraction_angle) if isinstance(sol, asin))


def calculate(incoming_intersection_y: float) -> tuple[Matrix, Matrix, Matrix] | None:
    """Calculates intersections with the lens and with the optical axis `y = 0`"""

    # Point of intersection of the incoming ray with the left side of the lens
    (incoming_intersection_x,) = (
        x for x in solve(left_side_eqn.subs(y, incoming_intersection_y), x) if x < 0)
    incoming_intersection_point = make_point(incoming_intersection_x, incoming_intersection_y)

    # Unit normal on the left surface of the lens at `incoming_intersection_point`
    incoming_unit_normal = normalize(incoming_intersection_point - left_side_center)

    incoming_ray_unit_vector = make_point(1, 0)

    sine_of_incoming_incidence_angle = sin_between_unit_vectors(
        incoming_unit_normal,
        incoming_ray_unit_vector,
    )

    incoming_refraction_angle = refraction_angle_expr.subs({
        snells_law.incidence_refractive_index: media_refractive_index,
        sin(snells_law.incidence_angle): sine_of_incoming_incidence_angle,
        snells_law.resulting_refractive_index: lens_refractive_index,
    })

    ray_unit_vector_in_lens = rotate(
        sign(incoming_intersection_y) * incoming_refraction_angle,
        -incoming_unit_normal,
    )
    ray_in_lens_eqn = Eq((p - incoming_intersection_point).cross(ray_unit_vector_in_lens)[2], 0)

    # Point of intersection of the ray in the lens and the right side of the lens.
    (outgoing_intersection_point,) = (make_point(sol[x], sol[y])
        for sol in solve((right_side_eqn, ray_in_lens_eqn), (x, y), dict=True)
        if sol[x] > 0)

    # Unit normal on the right surface of the lens at `outgoing_intersection_point`
    outgoing_unit_normal = normalize(outgoing_intersection_point - right_side_center)

    sine_of_outgoing_incidence_angle = sin_between_unit_vectors(
        outgoing_unit_normal,
        ray_unit_vector_in_lens,
    )

    outgoing_refraction_angle = refraction_angle_expr.subs({
        snells_law.incidence_refractive_index: lens_refractive_index,
        sin(snells_law.incidence_angle): sine_of_outgoing_incidence_angle,
        snells_law.resulting_refractive_index: media_refractive_index,
    })
    if im(outgoing_refraction_angle) != 0:
        # Case of total internal reflection
        return None

    outgoing_ray_unit_vector = rotate(
        -sign(incoming_intersection_y) * outgoing_refraction_angle,
        outgoing_unit_normal,
    )

    outgoing_ray_eqn = Eq(
        (p - outgoing_intersection_point).cross(outgoing_ray_unit_vector).dot(e_z),
        0,
    )

    # Point of intersection of the outgoing ray with the optical axis
    optical_axis_intersection_x = solve(outgoing_ray_eqn.subs(y, 0), x)[0]
    optical_axis_intersection_point = make_point(optical_axis_intersection_x, 0)

    return incoming_intersection_point, outgoing_intersection_point, optical_axis_intersection_point


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
        (effective_focus_length, y),
        (y, -maxval, maxval),
        line_color="pink",
        label="effective focus plane",
        show=False,
    )
    base_plot.extend(subplot)

    phi = sym_symbols("phi", real=True)
    left_max_angle: Expr = asin(lens_radius / left_radius)
    right_max_angle: Expr = asin(lens_radius / right_radius)

    eqn1_parametric = rotate(phi, make_point(left_radius, 0)) + left_side_center
    subplot = plot_parametric(
        eqn1_parametric[:-1],
        (phi, pi - left_max_angle, pi + left_max_angle),
        line_color="blue",
        label="",
        show=False,
    )
    base_plot.extend(subplot)

    eqn2_parametric = rotate(phi, make_point(right_radius, 0)) + right_side_center
    subplot = plot_parametric(
        eqn2_parametric[:-1],
        (phi, -right_max_angle, right_max_angle),
        line_color="blue",
        label="",
        show=False,
    )
    base_plot.extend(subplot)

    for p_in, p_out, p_f in data:
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
