from collections import namedtuple
from pytest import fixture, raises
from sympy import sin, cos, pi
from symplyphysics import units, Symbol, assert_equal, errors
from symplyphysics.core.vectors import vector_diff, VectorNorm, VectorCross
from symplyphysics.core.coordinate_systems import CARTESIAN, CoordinateVector
from symplyphysics.core.coordinate_systems.vector import as_coordinate_vector
from symplyphysics.core.coordinate_systems.point import AppliedPoint
from symplyphysics.core.coordinate_systems.curve import Curve
from symplyphysics.reorganized.mathematics.fields import flux_is_integral_along_flat_curve as law

Args = namedtuple("Args", "f n c v1 v2")

x, y, _ = CARTESIAN.base_scalars


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = 2 * units.meter
    b = 1 * units.meter

    v = Symbol("v", positive=True)

    pt = AppliedPoint([a * cos(v), b * sin(v), 0], CARTESIAN)
    c = Curve(v, pt)

    r = CoordinateVector(
        CARTESIAN.extract_position_vector_components(pt.coordinates.values()),
        CARTESIAN,
    )

    dr_dv = as_coordinate_vector(vector_diff(r, v))
    t = as_coordinate_vector(dr_dv / VectorNorm(dr_dv)).simplify()  # unit tangent vector
    e_z = CoordinateVector([0, 0, 1], CARTESIAN)
    n = as_coordinate_vector(VectorCross(t, e_z)).simplify()  # unit outward normal

    v1 = 0
    v2 = 2 * pi

    f = CoordinateVector(
        [1 * units.volt / units.meter, (x + y) * units.volt / units.meter**2, 0],
        CARTESIAN,
    )

    return Args(f=f, n=n, c=c, v1=v1, v2=v2)


def test_law(test_args: Args) -> None:
    result = law.calculate_flux(test_args.f, test_args.n, test_args.c, test_args.v1, test_args.v2)
    assert_equal(result, 2 * pi * units.volt, relative_tolerance=3e-3)


def test_bad_bounds(test_args: Args) -> None:
    bad_parameter = 2 * units.second

    with raises(errors.UnitsError):
        law.calculate_flux(test_args.f, test_args.n, test_args.c, bad_parameter, test_args.v2)
    with raises(errors.UnitsError):
        law.calculate_flux(test_args.f, test_args.n, test_args.c, test_args.v1, bad_parameter)
    with raises(errors.UnitsError):
        law.calculate_flux(test_args.f, test_args.n, test_args.c, bad_parameter, 2 * bad_parameter)
