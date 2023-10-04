from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import S, Expr, sin, cos, pi, sqrt
from symplyphysics import convert_to
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.field_point import FieldPoint
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.laws.fields import flux_is_integral_across_surface as flux_def


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    Args = namedtuple("Args", ["C"])
    return Args(C=C)


def test_basic_flux(test_args):
    field = VectorField(lambda point: [point.x, point.y], test_args.C)
    # flux over circle of radius = 2
    curve = [2 * cos(flux_def.parameter), 2 * sin(flux_def.parameter)]
    result = flux_def.calculate_flux(field, curve, (0, 2 * pi))
    assert convert_to(result, S.One).evalf(4) == approx((8 * pi).evalf(4), 0.001)


def test_ellipse_flux(test_args):
    field = VectorField(lambda point: [4 * point.x, point.y], test_args.C)
    # flux over ellipse with axes 2 and 1
    curve = [2 * cos(flux_def.parameter), sin(flux_def.parameter)]
    result = flux_def.calculate_flux(field, curve, (0, 2 * pi))
    assert convert_to(result, S.One).evalf(4) == approx((10 * pi).evalf(4), 0.001)


def test_three_dimensional_flux(test_args):
    field = VectorField(lambda point: [point.x, point.y], test_args.C)
    # trajectory is upwards helix
    helix = [
        cos(flux_def.parameter),
        sin(flux_def.parameter),
        flux_def.parameter
    ]
    with raises(ValueError):
        flux_def.calculate_flux(field, helix, (0, 2 * pi))


def _distance(point: FieldPoint) -> Expr:
    return sqrt(point.x**2 + point.y**2 + point.z**2)


def test_gravitational_field_flux(test_args):
    # gravitational field also has a common multiplier of -G*m1*m2.
    field = VectorField(
        lambda point: [
        point.x / _distance(point)**3, point.y / _distance(point)**3, point.z / _distance(point)**3], test_args.C)
    curve = [cos(flux_def.parameter), sin(flux_def.parameter)]
    result = flux_def.calculate_flux(field, curve, (0, 2 * pi))
    assert result > 0
    