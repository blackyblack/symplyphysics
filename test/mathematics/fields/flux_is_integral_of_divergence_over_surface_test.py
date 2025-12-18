from collections import namedtuple
from pytest import fixture, raises
from sympy import sin, cos, pi
from symplyphysics import units, Symbol, assert_equal, errors
from symplyphysics.core.coordinate_systems import CARTESIAN, CoordinateVector
from symplyphysics.core.coordinate_systems.point import AppliedPoint
from symplyphysics.core.coordinate_systems.surface import Surface
from symplyphysics.mathematics.fields import flux_is_integral_of_divergence_over_surface as law

Args = namedtuple("Args", "f s u1 u2 v1 v2")

x, y, _ = CARTESIAN.base_scalars


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = 2 * units.meter
    b = 1 * units.meter

    u = Symbol("v", positive=True)
    v = Symbol("v", positive=True)

    pt = AppliedPoint([a * u * cos(v), b * u * sin(v), 0], CARTESIAN)
    s = Surface((u, v), pt)

    u1 = 0
    u2 = 1
    v1 = 0
    v2 = 2 * pi

    f = CoordinateVector(
        [1 * units.volt / units.meter, (x + y) * units.volt / units.meter**2, 0],
        CARTESIAN,
    )

    return Args(f=f, s=s, u1=u1, u2=u2, v1=v1, v2=v2)


def test_law(test_args: Args) -> None:
    result = law.calculate_flux(test_args.f, test_args.s, test_args.u1, test_args.u2, test_args.v1,
        test_args.v2)
    assert_equal(result, 2 * pi * units.volt)


def test_bad_bounds(test_args: Args) -> None:
    bad_parameter = 1 * units.ampere

    with raises(errors.UnitsError):
        law.calculate_flux(test_args.f, test_args.s, bad_parameter, test_args.u2, test_args.v1,
            test_args.v2)
    with raises(errors.UnitsError):
        law.calculate_flux(test_args.f, test_args.s, test_args.u1, bad_parameter, test_args.v1,
            test_args.v2)
    with raises(errors.UnitsError):
        law.calculate_flux(test_args.f, test_args.s, bad_parameter, 2 * bad_parameter, test_args.v1,
            test_args.v2)
    with raises(errors.UnitsError):
        law.calculate_flux(test_args.f, test_args.s, test_args.u1, test_args.u2, bad_parameter,
            test_args.v2)
    with raises(errors.UnitsError):
        law.calculate_flux(test_args.f, test_args.s, test_args.u1, test_args.u2, test_args.v1,
            bad_parameter)
    with raises(errors.UnitsError):
        law.calculate_flux(test_args.f, test_args.s, test_args.u1, test_args.u2, bad_parameter,
            2 * bad_parameter)
