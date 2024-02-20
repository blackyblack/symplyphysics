from collections import namedtuple
from typing import Sequence
from pytest import fixture, mark, raises
from sympy import Expr, cos, pi, sin, sqrt, Symbol as SymSymbol, sympify
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.analysis import circulation_along_curve, circulation_along_surface_boundary, flux_across_curve, flux_across_surface, flux_across_surface_boundary, flux_across_volume_boundary
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.points.cartesian_point import CartesianPoint
from symplyphysics.core.points.cylinder_point import CylinderPoint
from symplyphysics.core.approx import approx_equal_numbers

Args = namedtuple("Args", ["C", "parameter1", "parameter2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C = CoordinateSystem()
    parameter1 = SymSymbol("parameter1")
    parameter2 = SymSymbol("parameter2")
    return Args(C=C, parameter1=parameter1, parameter2=parameter2)


def test_basic_circulation_along_curve(test_args: Args) -> None:
    field = VectorField(lambda point: [point.y, 0, point.x + point.z], test_args.C)
    curve = [cos(test_args.parameter1), sin(test_args.parameter1)]
    result = circulation_along_curve(field, curve, (test_args.parameter1, 0, pi / 2))
    result_sym = sympify(result)
    assert approx_equal_numbers(result_sym.evalf(4), (-pi / 4).evalf(4))


def test_circle_circulation(test_args: Args) -> None:
    field = VectorField(lambda point: [-point.y, point.x])
    radius = SymSymbol("radius")
    circle = [radius * cos(test_args.parameter1), radius * sin(test_args.parameter1)]
    result = circulation_along_curve(field, circle, (test_args.parameter1, 0, 2 * pi))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (2 * pi * radius**2).evalf(4))


def test_implicit_circulation_along_curve(test_args: Args) -> None:
    field = VectorField(lambda point: [point.y, -point.x, 0], test_args.C)
    # circle function is: x**2 + y**2 = 9
    # parametrize by circulation_def.parameter
    circle = [3 * cos(test_args.parameter1), 3 * sin(test_args.parameter1)]
    result = circulation_along_curve(field, circle, (test_args.parameter1, 0, 2 * pi))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (-18 * pi).evalf(4))
    # now try to define trajectory without parametrization
    # parametrized solution uses angle [0, 2*pi] that corresponds to the counter-clockwise direction
    # so we should integrate in the same direction: [r, -r] for upper part of the circle and [-r, r] for lower
    # y = sqrt(9 - x**2) for upper part of the circle
    # y = -sqrt(9 - x**2) for lower part of the circle
    circle_implicit_up = [test_args.parameter1, sqrt(9 - test_args.parameter1**2)]
    result_up = circulation_along_curve(field, circle_implicit_up, (test_args.parameter1, 3, -3))
    circle_implicit_down = [test_args.parameter1, -sqrt(9 - test_args.parameter1**2)]
    result_down = circulation_along_curve(field, circle_implicit_down,
        (test_args.parameter1, -3, 3))
    result_up_sym = sympify(result_up)
    result_down_sym = sympify(result_down)
    approx_equal_numbers(result_up_sym.evalf(4) + result_down_sym.evalf(4), (-18 * pi).evalf(4))


def test_orthogonal_movement_circulation_along_curve(test_args: Args) -> None:
    field = VectorField(lambda point: [point.y, -point.x, 1], test_args.C)
    # trajectory is upwards helix
    helix = [cos(test_args.parameter1), sin(test_args.parameter1), test_args.parameter1]
    result = circulation_along_curve(field, helix, (test_args.parameter1, 0, 2 * pi))
    assert result == 0
    # trajectory is upwards straight line
    trajectory_vertical = [1, 0, test_args.parameter1]
    result = circulation_along_curve(field, trajectory_vertical, (test_args.parameter1, 0, 2 * pi))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (2 * pi).evalf(4))


def test_basic_circulation_along_surface_boundary(test_args: Args) -> None:
    field = VectorField(lambda point: [point.y, 0, point.x + point.z], test_args.C)
    # cylinder
    surface = [
        test_args.parameter1 * cos(test_args.parameter2),
        test_args.parameter1 * sin(test_args.parameter2)
    ]
    result = circulation_along_surface_boundary(field, surface, (test_args.parameter1, 0, 1),
        (test_args.parameter2, 0, pi / 2))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (-pi / 4).evalf(4))
    # circle, which is a cylinder's boundary
    curve = [cos(test_args.parameter1), sin(test_args.parameter1)]
    result_from_boundary = circulation_along_curve(field, curve, (test_args.parameter1, 0, pi / 2))
    # verify that 'circulation_along_curve()' result is the same as when using
    # 'circulation_along_surface_boundary()'
    result_from_boundary_sym = sympify(result_from_boundary)
    approx_equal_numbers(result_from_boundary_sym.evalf(4), result_sym.evalf(4))


def test_cone_circulation_along_surface_boundary(test_args: Args) -> None:
    field = VectorField(lambda point: [point.y, -point.x, 0], test_args.C)
    # circle function is: x**2 + y**2 = 9
    # from circulation_is_integral_along_curve_test we got circulation -18 * pi
    # let's check with cone surface
    cone = [
        3 * test_args.parameter1 * cos(test_args.parameter2),
        3 * test_args.parameter1 * sin(test_args.parameter2), test_args.parameter1
    ]
    result = circulation_along_surface_boundary(field, cone, (test_args.parameter1, 0, 1),
        (test_args.parameter2, 0, 2 * pi))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (-18 * pi).evalf(4))


# non-surface trajectories (eg circle) result in zero circulation
def test_circle_circulation_along_surface_boundary(test_args: Args) -> None:
    field = VectorField(lambda point: [point.y, point.x], test_args.C)
    surface = [cos(test_args.parameter1), sin(test_args.parameter1)]
    result = circulation_along_surface_boundary(field, surface, (test_args.parameter1, 0, 1),
        (test_args.parameter2, 0, pi / 2))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), 0)


def test_basic_flux_across_curve(test_args: Args) -> None:
    field = VectorField(lambda point: [point.x, point.y], test_args.C)
    # flux over circle of radius = 2
    curve = [2 * cos(test_args.parameter1), 2 * sin(test_args.parameter1)]
    result = flux_across_curve(field, curve, (test_args.parameter1, 0, 2 * pi))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (8 * pi).evalf(4))


def test_ellipse_flux(test_args: Args) -> None:
    field = VectorField(lambda point: [4 * point.x, point.y], test_args.C)
    # flux over ellipse with axes 2 and 1
    curve = [2 * cos(test_args.parameter1), sin(test_args.parameter1)]
    result = flux_across_curve(field, curve, (test_args.parameter1, 0, 2 * pi))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (10 * pi).evalf(4))


# 3-dimensional curves are not supported
def test_helix_flux_across_curve(test_args: Args) -> None:
    field = VectorField(lambda point: [point.x, point.y], test_args.C)
    # trajectory is upwards helix
    helix = [cos(test_args.parameter1), sin(test_args.parameter1), test_args.parameter1]
    with raises(ValueError):
        flux_across_curve(field, helix, (test_args.parameter1, 0, 2 * pi))


def test_basic_flux_across_surface(test_args: Args) -> None:
    field = VectorField(lambda point: [point.x, point.y, 0], test_args.C)
    # flux over hyperboloid, between planes z = -2 and z = 1
    curve = [
        sqrt(test_args.parameter1**2 + 1) * cos(test_args.parameter2),
        sqrt(test_args.parameter1**2 + 1) * sin(test_args.parameter2), -test_args.parameter1
    ]
    result = flux_across_surface(field, curve, (test_args.parameter1, -2, 1),
        (test_args.parameter2, 0, 2 * pi))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (12 * pi).evalf(4))


def _distance(point: CartesianPoint) -> Expr:
    return sqrt(point.x**2 + point.y**2 + point.z**2)


@mark.skip("takes too much time to calculate the result")
def test_gravitational_field_flux(test_args: Args) -> None:
    # gravitational field also has a common multiplier of -G * M
    field = VectorField(
        lambda point: [
        point.x / _distance(point)**2, point.y / _distance(point)**2, point.z / _distance(point)**2
        ], test_args.C)
    trajectory = [
        cos(test_args.parameter1) * sin(test_args.parameter2),
        sin(test_args.parameter1) * sin(test_args.parameter2),
        cos(test_args.parameter2)
    ]
    result = flux_across_surface(field, trajectory, (test_args.parameter1, 0, 2 * pi),
        (test_args.parameter2, 0, pi))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (4 * pi).evalf(4))
    # if -G * M multiplier is applied, we will get Gauss law for gravity:
    # https://en.wikipedia.org/wiki/Gauss%27s_law_for_gravity


# it works when field divergence is constant (no free variables), or when surface is not
# parametrized
def test_basic_flux_across_surface_boundary(test_args: Args) -> None:
    field = VectorField(lambda point: [point.x - point.y, point.x + point.y], test_args.C)
    # flux over circle of radius 2
    circle = [
        test_args.parameter1 * cos(test_args.parameter2),
        test_args.parameter1 * sin(test_args.parameter2)
    ]
    result = flux_across_surface_boundary(field, circle, (test_args.parameter1, 0, 2),
        (test_args.parameter2, 0, 2 * pi))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (8 * pi).evalf(4))


# non parametrized surface works for any divergence, but is harder to manipulate with
def test_non_parametrized_flux_across_surface_boundary(test_args: Args) -> None:
    field = VectorField(lambda point: [point.x**2, point.y], test_args.C)
    x = field.coordinate_system.coord_system.base_scalars()[0]
    y = field.coordinate_system.coord_system.base_scalars()[1]
    circle_implicit = [x, y]
    # flux over circle of radius 3
    result = flux_across_surface_boundary(field, circle_implicit,
        (x, -sqrt(9 - y**2), sqrt(9 - y**2)), (y, -3, 3))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (9 * pi).evalf(4))


def test_basic_flux_across_volume_boundary(test_args: Args) -> None:
    field = VectorField(lambda point: [point.x**2 / 2, point.y * point.z, -point.x * point.z],
        test_args.C)
    x = field.coordinate_system.coord_system.base_scalars()[0]
    y = field.coordinate_system.coord_system.base_scalars()[1]
    result = flux_across_volume_boundary(field, (-2, 2), (-sqrt(4 - x**2), sqrt(4 - x**2)),
        (0, sqrt(4 - x**2 - y**2)))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (4 * pi).evalf(4))


# A vector field exists in the region between two concentric cylindrical surfaces defined by
# r = 1 and r = 2, with both cylinders extending between z = 0 and z = 5
def test_basic_flux_across_sphere_boundary() -> None:
    B = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)

    def field_function(p: CylinderPoint) -> Sequence[ScalarValue]:
        return [p.radius**3, 0, 0]

    field = VectorField(field_function, B)
    result = flux_across_volume_boundary(field, (1, 2), (0, 2 * pi), (0, 5))
    result_sym = sympify(result)
    approx_equal_numbers(result_sym.evalf(4), (150 * pi).evalf(4))
