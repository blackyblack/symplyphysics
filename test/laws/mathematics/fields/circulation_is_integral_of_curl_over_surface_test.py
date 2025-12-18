from collections import namedtuple
from pytest import fixture, raises
from sympy import sin, sqrt, pi, Symbol as SymSymbol
from symplyphysics import assert_equal, units
from symplyphysics.core.errors import UnitsError
from symplyphysics.core.coordinate_systems import CYLINDRICAL, CoordinateVector
from symplyphysics.core.coordinate_systems.point import AppliedPoint
from symplyphysics.core.coordinate_systems.surface import Surface
from symplyphysics.laws.mathematics.fields import circulation_is_integral_of_curl_over_surface as law

rho, phi, z = CYLINDRICAL.base_scalars

Args = namedtuple("Args", "f s u1 u2 v1 v2")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    u = SymSymbol("u", positive=True)
    v = SymSymbol("v", positive=True)

    pt = AppliedPoint([u * units.meter, v, u * sin(v) * units.meter], CYLINDRICAL)
    s = Surface((u, v), pt)

    f = CoordinateVector(
        [z * units.volt / units.meter**2, 0,
        sin(phi / 5) * units.volt / units.meter],
        CYLINDRICAL,
        pt,
    )

    u1 = 0
    u2 = 1
    v1 = 0
    v2 = 2 * pi

    return Args(f=f, s=s, u1=u1, u2=u2, v1=v1, v2=v2)


def test_law(test_args: Args) -> None:
    result = law.calculate_circulation(test_args.f, test_args.s, test_args.u1, test_args.u2,
        test_args.v1, test_args.v2)
    assert_equal(result, 5 * (sqrt(5) - 5) / 96 * units.volt)


def test_bad_bounds(test_args: Args) -> None:
    bad_parameter = 1 * units.ampere

    with raises(UnitsError):
        law.calculate_circulation(test_args.f, test_args.s, bad_parameter, test_args.u2,
            test_args.v1, test_args.v2)
    with raises(UnitsError):
        law.calculate_circulation(test_args.f, test_args.s, test_args.u1, bad_parameter,
            test_args.v1, test_args.v2)
    with raises(UnitsError):
        law.calculate_circulation(test_args.f, test_args.s, bad_parameter, bad_parameter,
            test_args.v1, test_args.v2)
    with raises(UnitsError):
        law.calculate_circulation(test_args.f, test_args.s, test_args.u1, test_args.u2,
            bad_parameter, test_args.v2)
    with raises(UnitsError):
        law.calculate_circulation(test_args.f, test_args.s, test_args.u1, test_args.u2,
            test_args.v1, bad_parameter)
    with raises(UnitsError):
        law.calculate_circulation(test_args.f, test_args.s, test_args.u1, test_args.u2,
            bad_parameter, bad_parameter)
