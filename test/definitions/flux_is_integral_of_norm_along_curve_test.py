from collections import namedtuple
from pytest import approx, fixture
from sympy import S, sin, cos, pi
from symplyphysics import (
    convert_to,)
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.definitions import flux_is_integral_of_norm_along_curve as flux_def


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
    assert float(convert_to(result, S.One).evalf(4)) == approx(float((8 * pi).evalf(4)), 0.001)
