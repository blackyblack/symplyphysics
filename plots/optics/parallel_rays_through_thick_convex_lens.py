#!/usr/bin/env python3
"""
Plot the passage of rays parallel to the optical axis through a thick convex lens.
"""

from __future__ import annotations
from typing import Any
from math import isnan, isinf
from argparse import ArgumentParser
import sys
from sympy import (ImmutableMatrix, symbols as sym_symbols, Eq, solve, Expr, sqrt, rot_axis3, sin,
    pi, asin, Point2D, sign, im, plot, plot_parametric)
from symplyphysics.core.coordinate_systems import CoordinateVector, CARTESIAN
from symplyphysics.core.coordinate_systems.vector import as_coordinate_vector, component
from symplyphysics.core.vectors import VectorNorm, VectorCross, VectorDot
from symplyphysics.laws.optics import refraction_angle_from_environments as snells_law
from symplyphysics.core.geometry.line import two_point_function

# PART 1. Parsing command line arguments


def make_parser() -> ArgumentParser:
    parser = ArgumentParser()
    parser.add_argument(
        "-r1",
        action="store",
        default=2.0,
        dest="left_side_radius",
        help="radius of curvature of the left side of the lens",
    )
    parser.add_argument(
        "-r2",
        action="store",
        default=2.0,
        dest="right_side_radius",
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
left_side_radius = positive_finite_float(args.left_side_radius)

# Radius of curvature of the **right** side of the lens
right_side_radius = positive_finite_float(args.right_side_radius)

# Refractive index of the medium surrounding the lens
media_refractive_index = positive_finite_float(args.media_refractive_index)

# Refractive index of the material of the lens
lens_refractive_index = positive_finite_float(args.lens_refractive_index)

if media_refractive_index >= lens_refractive_index:
    MESSAGE = "Expect the refraction index of the surrounding medium to be smaller than that of the lens."
    raise ValueError(MESSAGE)

# PART 2. Calculations


def make_point(x_: Any, y_: Any) -> CoordinateVector:
    """Embeds a 2-dimensional point `(x, y)` into the 3-dimensional plane `z = 0`."""

    return CoordinateVector([x_, y_, 0], CARTESIAN)


# NOTE Refer to the figure at `../../img/parallel_rays_through_thick_convex_lens.svg`

# Radius `GH` of the lens, i.e. the radius at the base of the [spherical caps](https://en.wikipedia.org/wiki/Spherical_cap) which compose the lens.
lens_radius = 1  # pylint: disable=invalid-name

# `C1C2 = C1H + C2H`
# Pythagoras theorem: `C1H^2 + GH^2 = R1^2`
# Pythagoras theorem: `C2H^2 + GH^2 = R2^2`
distance_between_centers = (sqrt(left_side_radius**2 - lens_radius**2) +
    sqrt(right_side_radius**2 - lens_radius**2))

# `C2H` using Pythagoras theorem in the triangle `C2HG`
left_side_center_to_origin = sqrt(left_side_radius**2 - lens_radius**2)

# `C1H` using Pythagoras theorem in the triangle `C1HG`
right_side_center_to_origin = sqrt(right_side_radius**2 - lens_radius**2)

# `C1C2 = C1H + C2H` since `C1, H, C2` lie on a straight line.
distance_between_centers = left_side_center_to_origin + right_side_center_to_origin

# `C2E2 = C1C2 - C1E2`
left_side_center_to_lens = distance_between_centers - right_side_radius

# `C1E1 = C1C2 - C2E1`
right_side_center_to_lens = distance_between_centers - left_side_radius

# Thickness `E1E2` of the lens: `E1E2 = C1C2 - C1E1 - C2E2`
lens_thickness = distance_between_centers - left_side_center_to_lens - right_side_center_to_lens

# Position vector of the center of the circle making up the **left** side of the lens
left_side_center = make_point(+left_side_center_to_origin, 0)

# Position vector of the center of the circle making up the **right** side of the lens
right_side_center = make_point(-right_side_center_to_origin, 0)

# Represents a point on the plane. To be used in equations as free variables.
x, y = sym_symbols("x y", real=True)
p = make_point(x, y)


def sq_norm(v: Expr) -> Expr:
    """Returns the square of the norm of `v`."""
    return VectorNorm(as_coordinate_vector(v)).doit()**2


def circle_equation(center: CoordinateVector, radius: Expr) -> Eq:
    # A circle is a set of points positioned at an equal distance (`radius`) from a given point (`center`).
    # In Cartesian coordinates, this equation is `(x - x_center)^2 + (y - y_center)^2 = radius^2`
    return Eq(sq_norm(p - center), radius**2)


# Equation of the **left** side of the lens
left_side_eqn = circle_equation(left_side_center, left_side_radius)

# Equation of the **right** side of the lens
right_side_eqn = circle_equation(right_side_center, right_side_radius)
# right_side_eqn = Eq(sq_norm(p - right_side_center), right_side_radius**2)


def normalize(v: Expr) -> CoordinateVector:
    """Returns the unit direction vector co-directed with the given vector."""
    return as_coordinate_vector(v / sqrt(sq_norm(v)))


e_z = CoordinateVector([0, 0, 1], CARTESIAN)


def sin_between_unit_vectors(u1: CoordinateVector, u2: CoordinateVector) -> Expr:
    """Assumes unit input vectors with a zero third component."""
    return abs(VectorDot(e_z, VectorCross(u1, u2))).doit()


def rotate(theta: Expr, v: Expr) -> CoordinateVector:
    """
    Rotates vector `v` anti-clockwise at angle `theta` if `theta > 0`, and clockwise otherwise.

    >>> assert rotate(pi / 2, make_point(1, 0)) == make_point(0, 1)
    >>> assert rotate(-pi / 2, make_point(1, 1)) == make_point(1, -1)
    """

    v = as_coordinate_vector(v)

    old_components = ImmutableMatrix([component(v, 0), component(v, 1), component(v, 2)])
    new_components = rot_axis3(-theta) * old_components

    return CoordinateVector(new_components, v.system, v.point)


def line_eqn(point_through: Expr, direction_vector: Expr) -> Expr:
    # If `p` is a free position vector, `q` is the position vector of a point on the line, and `u`
    # is the unit direction vector of the line, then `p - q` would be parallel to `u` if and only
    # if `p` lies on the line as well. Therefore, the cross product between `p - q` and `u` would
    # be `0`. Since all of `p, q, u` lie in the `xy`-plane, the cross product would be a vector
    # parallel to the `z`-axis. Here, `q` is `point_through` and `u` is `direction_vector`.

    return VectorDot(e_z, VectorCross(direction_vector, p - point_through)).doit()


# The other solution is `pi - sol`; we don't need it since the refraction angle cannot exceed `pi / 2`.
(refraction_angle_expr,) = (
    sol for sol in solve(snells_law.law, snells_law.refraction_angle) if isinstance(sol, asin))


def calculate(
    incoming_intersection_y: float,
) -> tuple[CoordinateVector, CoordinateVector, CoordinateVector] | None:
    """Calculates intersections with the lens and the optical axis `y = 0`"""

    # Point `B1` of intersection of the incoming ray with the left side of the lens
    (incoming_intersection_x,) = (
        x for x in solve(left_side_eqn.subs(y, incoming_intersection_y), x) if x < 0)
    incoming_intersection_point = make_point(incoming_intersection_x, incoming_intersection_y)

    # Unit normal `B1N1` on the left surface of the lens at `incoming_intersection_point`
    incoming_unit_normal = normalize(incoming_intersection_point - left_side_center)

    # Unit vector parallel to `A1B1`
    incoming_ray_unit_vector = make_point(1, 0)

    # Sine of angle `N1B1A1`
    sine_of_incoming_incidence_angle = sin_between_unit_vectors(incoming_unit_normal,
        incoming_ray_unit_vector)

    # Angle `M1B1B2`
    incoming_refraction_angle = refraction_angle_expr.subs({
        snells_law.incidence_refractive_index: media_refractive_index,
        sin(snells_law.incidence_angle): sine_of_incoming_incidence_angle,
        snells_law.resulting_refractive_index: lens_refractive_index,
    })

    # Unit vector parallel to `B1B2`. We can obtain it by rotating `B1M1` by `incoming_refraction_angle`.
    # The rotating is counterclockwise in the "top" part of the lens, else clockwise (see figure).
    ray_unit_vector_in_lens = rotate(
        sign(incoming_intersection_y) * incoming_refraction_angle,
        -incoming_unit_normal,
    )
    ray_in_lens_eqn = line_eqn(incoming_intersection_point, ray_unit_vector_in_lens)

    # Point `B2` of intersection of the ray in the lens and the right side of the lens.
    (outgoing_intersection_point,) = (make_point(sol[x], sol[y])
        for sol in solve((right_side_eqn, ray_in_lens_eqn), (x, y), dict=True)
        if sol[x] > 0)

    # Unit normal `B2N2` on the right surface of the lens at `outgoing_intersection_point`
    outgoing_unit_normal = normalize(outgoing_intersection_point - right_side_center)

    # Sine of angle `B1B2M2`
    sine_of_outgoing_incidence_angle = sin_between_unit_vectors(outgoing_unit_normal,
        ray_unit_vector_in_lens)

    # Angle `N2B2A2`
    outgoing_refraction_angle = refraction_angle_expr.subs({
        snells_law.incidence_refractive_index: lens_refractive_index,
        sin(snells_law.incidence_angle): sine_of_outgoing_incidence_angle,
        snells_law.resulting_refractive_index: media_refractive_index,
    })

    # Total internal reflection: there is no outgoing ray.
    if im(outgoing_refraction_angle) != 0:
        return None

    # Unit vector parallel to `B2A2`. We can obtain it by rotating `B2N2` by angle `outgoing_refraction_angle`.
    # The rotation is clockwise in the "top" part of the lens, else counterclockwise (see figure).
    outgoing_ray_unit_vector = rotate(
        -sign(incoming_intersection_y) * outgoing_refraction_angle,
        outgoing_unit_normal,
    )

    # Equation of the line `B2A2`
    outgoing_ray_eqn = line_eqn(outgoing_intersection_point, outgoing_ray_unit_vector)

    # Point `A2` of intersection of the outgoing ray with the optical axis
    optical_axis_intersection_x = solve(outgoing_ray_eqn.subs(y, 0), x)[0]
    optical_axis_intersection_point = make_point(optical_axis_intersection_x, 0)

    # Return `B1, B2, A2`
    return (incoming_intersection_point, outgoing_intersection_point,
        optical_axis_intersection_point)


# PART 3. Plotting


def to_point2d(v: Expr) -> Point2D:
    v = as_coordinate_vector(v)

    v0 = component(v, 0)
    v1 = component(v, 1)

    return Point2D(v0, v1)


def make_line_eqn(v1: Expr, v2: Expr) -> Expr:
    return two_point_function(to_point2d(v1), to_point2d(v2), x)


def display() -> None:
    y_ins = (-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9)

    data: list[tuple[CoordinateVector, CoordinateVector, CoordinateVector]] = []
    for y_in in y_ins:
        t = calculate(y_in)
        if not t:
            continue
        data.append(t)

    maxval = max(max((component(t[2], 0) for t in data), default=1), 1) * 1.1

    base_plot = plot(
        title="Parallel rays passing through a thick convex lens",
        xlabel="",
        ylabel="",
        ylim=(-maxval, maxval),
        xlim=(-maxval, maxval),
        size=(8, 8),
        show=False,
    )

    phi = sym_symbols("phi", real=True)

    # Angle `GC2C1` (see figure)
    left_max_angle: Expr = asin(lens_radius / left_side_radius)

    # Angle `GC1C2` (see figure)
    right_max_angle: Expr = asin(lens_radius / right_side_radius)

    # Rotate the radius-vector relative to the circle center and then shift the circle center.
    eqn1_parametric = as_coordinate_vector(
        rotate(phi, make_point(left_side_radius, 0)) + left_side_center)
    subplot = plot_parametric(
        eqn1_parametric.components[:-1],
        (phi, pi - left_max_angle, pi + left_max_angle),
        line_color="blue",
        label="",
        show=False,
    )
    base_plot.extend(subplot)

    # Rotate the radius-vector relative to the circle center and then shift the circle center.
    eqn2_parametric = as_coordinate_vector(
        rotate(phi, make_point(right_side_radius, 0)) + right_side_center)
    subplot = plot_parametric(
        eqn2_parametric.components[:-1],
        (phi, -right_max_angle, right_max_angle),
        line_color="blue",
        label="",
        show=False,
    )
    base_plot.extend(subplot)

    for p_in, p_out, p_f in data:
        p_left = make_point(-maxval, p_in.components[1])
        subplot = plot(
            make_line_eqn(p_left, p_in),
            (x, -maxval, p_in.components[0]),
            line_color="green",
            label="",
            show=False,
        )
        base_plot.extend(subplot)

        subplot = plot(
            make_line_eqn(p_in, p_out),
            (x, p_in.components[0], p_out.components[0]),
            line_color="green",
            label="",
            show=False,
        )
        base_plot.extend(subplot)

        subplot = plot(
            make_line_eqn(p_out, p_f),
            (x, p_out.components[0], p_f.components[0]),
            line_color="green",
            label="",
            show=False,
        )
        base_plot.extend(subplot)

    base_plot.show()


display()
