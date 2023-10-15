from collections import namedtuple
from pytest import approx, fixture, mark, raises
from sympy import Expr, cos, pi, sin, sqrt, Symbol as SymSymbol
from sympy.vector import VectorZero
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.analysis import circulation_along_curve, circulation_over_surface, flux_across_curve, flux_across_surface
from symplyphysics.core.fields.field_point import FieldPoint
from symplyphysics.core.fields.operators import curl_operator
from symplyphysics.core.fields.vector_field import VectorField


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    parameter1 = SymSymbol("parameter1")
    parameter2 = SymSymbol("parameter2")
    Args = namedtuple("Args", ["C", "parameter1", "parameter2"])
    return Args(C=C, parameter1=parameter1, parameter2=parameter2)


def _distance(point: FieldPoint) -> Expr:
    return sqrt(point.x**2 + point.y**2 + point.z**2)


def test_basic_circulation_along_curve(test_args):
    field = VectorField(lambda point: [point.y, 0, point.x + point.z], test_args.C)
    curve = [cos(test_args.parameter1), sin(test_args.parameter1)]
    result = circulation_along_curve(field, curve, (test_args.parameter1, 0, pi / 2))
    assert result.evalf(4) == approx((-pi / 4).evalf(4), 0.001)


def test_circle_circulation(test_args):
    field = VectorField(lambda point: [-point.y, point.x])
    radius = SymSymbol("radius")
    circle = [radius * cos(test_args.parameter1), radius * sin(test_args.parameter1)]
    result = circulation_along_curve(field, circle, (test_args.parameter1, 0, 2 * pi))
    assert result.evalf(4) == approx((2 * pi * radius**2).evalf(4), 0.001)


def test_implicit_circulation_along_curve(test_args):
    field = VectorField(lambda point: [point.y, -point.x, 0], test_args.C)
    # circle function is: x**2 + y**2 = 9
    # parametrize by circulation_def.parameter
    circle = [3 * cos(test_args.parameter1), 3 * sin(test_args.parameter1)]
    result = circulation_along_curve(field, circle, (test_args.parameter1, 0, 2 * pi))
    assert result.evalf(4) == approx((-18 * pi).evalf(4), 0.001)
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
    assert (result_up.evalf(4) + result_down.evalf(4)) == approx((-18 * pi).evalf(4), 0.001)


def test_orthogonal_movement_circulation_along_curve(test_args):
    field = VectorField(lambda point: [point.y, -point.x, 1], test_args.C)
    # trajectory is upwards helix
    helix = [cos(test_args.parameter1), sin(test_args.parameter1), test_args.parameter1]
    result = circulation_along_curve(field, helix, (test_args.parameter1, 0, 2 * pi))
    assert result == 0
    # trajectory is upwards straight line
    trajectory_vertical = [1, 0, test_args.parameter1]
    result = circulation_along_curve(field, trajectory_vertical, (test_args.parameter1, 0, 2 * pi))
    assert result.evalf(4) == approx((2 * pi).evalf(4), 0.001)


def test_basic_circulation_over_surface(test_args):
    field = VectorField(lambda point: [point.y, 0, point.x + point.z], test_args.C)
    surface = [
        test_args.parameter1 * cos(test_args.parameter2),
        test_args.parameter1 * sin(test_args.parameter2)
    ]
    result = circulation_over_surface(field, surface, (test_args.parameter1, 0, 1),
        (test_args.parameter2, 0, pi / 2))
    assert result.evalf(4) == approx((-pi / 4).evalf(4), 0.001)


def test_cone_circulation(test_args):
    field = VectorField(lambda point: [point.y, -point.x, 0], test_args.C)
    # circle function is: x**2 + y**2 = 9
    # from circulation_is_integral_along_curve_test we got circulation -18 * pi
    # let's check with cone surface
    cone = [
        3 * test_args.parameter1 * cos(test_args.parameter2),
        3 * test_args.parameter1 * sin(test_args.parameter2), test_args.parameter1
    ]
    result = circulation_over_surface(field, cone, (test_args.parameter1, 0, 1),
        (test_args.parameter2, 0, 2 * pi))
    assert result.evalf(4) == approx((-18 * pi).evalf(4), 0.001)


# non-surface trajectories (eg circle) are not supported
def test_circle_circulation_over_surface(test_args):
    field = VectorField(lambda point: [point.y, point.x], test_args.C)
    surface = [cos(test_args.parameter1), sin(test_args.parameter1)]
    with raises(ValueError):
        circulation_over_surface(field, surface, (test_args.parameter1, 0, 1),
            (test_args.parameter2, 0, pi / 2))


def test_gravitational_field_is_conservative(test_args):
    # gravitational field also has a common multiplier of -G*M. It does not
    # affect conservative property of a field.
    field = VectorField(
        lambda point: [
        point.x / _distance(point)**3, point.y / _distance(point)**3, point.z / _distance(point)**3
        ], test_args.C)
    field_rotor = curl_operator(field)
    field_rotor_applied = field_rotor.apply_to_basis().to_sympy_vector()
    assert field_rotor_applied == VectorZero.zero


def test_basic_flux_across_curve(test_args):
    field = VectorField(lambda point: [point.x, point.y], test_args.C)
    # flux over circle of radius = 2
    curve = [2 * cos(test_args.parameter1), 2 * sin(test_args.parameter1)]
    result = flux_across_curve(field, curve, (test_args.parameter1, 0, 2 * pi))
    assert result.evalf(4) == approx((8 * pi).evalf(4), 0.001)


def test_ellipse_flux(test_args):
    field = VectorField(lambda point: [4 * point.x, point.y], test_args.C)
    # flux over ellipse with axes 2 and 1
    curve = [2 * cos(test_args.parameter1), sin(test_args.parameter1)]
    result = flux_across_curve(field, curve, (test_args.parameter1, 0, 2 * pi))
    assert result.evalf(4) == approx((10 * pi).evalf(4), 0.001)


# 3-dimensional curves are not supported
def test_helix_flux_across_curve(test_args):
    field = VectorField(lambda point: [point.x, point.y], test_args.C)
    # trajectory is upwards helix
    helix = [cos(test_args.parameter1), sin(test_args.parameter1), test_args.parameter1]
    with raises(ValueError):
        flux_across_curve(field, helix, (test_args.parameter1, 0, 2 * pi))


def test_basic_flux_across_surface(test_args):
    field = VectorField(lambda point: [point.x, point.y, 0], test_args.C)
    # flux over hyperboloid, between planes z = -2 and z = 1
    curve = [
        sqrt(test_args.parameter1**2 + 1) * cos(test_args.parameter2),
        sqrt(test_args.parameter1**2 + 1) * sin(test_args.parameter2), -test_args.parameter1
    ]
    result = flux_across_surface(field, curve, (test_args.parameter1, -2, 1),
        (test_args.parameter2, 0, 2 * pi))
    assert result.evalf(4) == approx((12 * pi).evalf(4), 0.001)


@mark.skip("takes too much time to calculate the result")
def test_gravitational_field_flux(test_args):
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
    assert result.evalf(4) == approx((4 * pi).evalf(4), 0.001)
    # if -G * M multiplier is applied, we will get Gauss law for gravity:
    # https://en.wikipedia.org/wiki/Gauss%27s_law_for_gravity
