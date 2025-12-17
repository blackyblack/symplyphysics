from collections import namedtuple
from pytest import fixture, raises
from sympy import sin, pi, sqrt
from symplyphysics import units, assert_equal, errors
from symplyphysics.core.coordinate_systems import CYLINDRICAL, CoordinateVector
from symplyphysics.core.coordinate_systems.point import AppliedPoint
from symplyphysics.core.coordinate_systems.curve import Curve
from symplyphysics.reorganized.mathematics.fields import circulation_is_integral_along_curve as law

Args = namedtuple("Args", "f c u1 u2")

rho, phi, z = CYLINDRICAL.base_scalars


@fixture(name="test_args")
def test_args_fixture() -> Args:
    pt = AppliedPoint([units.meter, phi, sin(phi) * units.meter], CYLINDRICAL)
    c = Curve(phi, pt)

    f = CoordinateVector(
        [z * units.volt / units.meter**2, 0,
        sin(phi / 5) * units.volt / units.meter],
        CYLINDRICAL,
        pt,
    )

    u1 = 0
    u2 = 2 * pi

    return Args(f=f, c=c, u1=u1, u2=u2)


def test_law(test_args: Args) -> None:
    result = law.calculate_circulation(test_args.f, test_args.c, test_args.u1, test_args.u2)
    assert_equal(result, 5 * (sqrt(5) - 5) / 96 * units.volt)


def test_bad_bounds(test_args: Args) -> None:
    bad_parameter = 1 * units.meter

    with raises(errors.UnitsError):
        law.calculate_circulation(test_args.f, test_args.c, 1, bad_parameter)
    with raises(errors.UnitsError):
        law.calculate_circulation(test_args.f, test_args.c, bad_parameter, 1)
    with raises(errors.UnitsError):
        law.calculate_circulation(test_args.f, test_args.c, bad_parameter, bad_parameter)
