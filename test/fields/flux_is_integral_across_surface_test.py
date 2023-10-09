from collections import namedtuple
from pytest import approx, fixture
from sympy import (S, Expr, sin, cos, pi, sqrt)
from symplyphysics import (
    SI,
    Quantity,
    convert_to,
    units,
)
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.field_point import FieldPoint
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.laws.fields import flux_is_integral_across_surface as flux_def


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    force_unit = Quantity(1 * units.newton)
    radius_unit = Quantity(1 * units.meter)
    Args = namedtuple("Args", ["C", "force_unit", "radius_unit"])
    return Args(C=C, force_unit=force_unit, radius_unit=radius_unit)


def test_basic_flux(test_args):
    field = VectorField(lambda point: [point.x, point.y, 0], test_args.C)
    # flux over hyperboloid, between planes z = -2 and z = 1
    curve = [sqrt(flux_def.parameter1**2 + 1) * cos(flux_def.parameter2), sqrt(flux_def.parameter1**2 + 1) * sin(flux_def.parameter2), -flux_def.parameter1]
    result = flux_def.calculate_flux(field, curve, (-2, 1), (0, 2 * pi))
    assert convert_to(result, S.One).evalf(4) == approx((12 * pi).evalf(4), 0.001)


def _distance(point: FieldPoint) -> Expr:
    return sqrt(point.x**2 + point.y**2 + point.z**2)


def test_gravitational_field_flux(test_args):
    field = VectorField(
        lambda point: [
        point.x / _distance(point)**2, point.y / _distance(point)**2, point.z / _distance(point)**2
        ], test_args.C)
    trajectory = [cos(flux_def.parameter1) * sin(flux_def.parameter2), sin(flux_def.parameter1) * sin(flux_def.parameter2), cos(flux_def.parameter2)]
    result = flux_def.calculate_flux(field, trajectory,
        (0, 2 * pi), (0, pi))
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
