from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity import semiminor_axis_of_elliptical_orbit_via_orbit_parameters as law

Args = namedtuple("Args", "s a m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    s = Quantity(2.23e9 * units.kilometer**2 / units.second)
    a = Quantity(units.astronomical_unit)
    m = Quantity(2.00e30 * units.kilogram)
    return Args(s=s, a=a, m=m)


def test_law(test_args: Args) -> None:
    result = law.calculate_semiminor_axis(test_args.s, test_args.a, test_args.m)
    assert_equal(result, 0.9999 * units.astronomical_unit, tolerance=2e-3)


def test_bad_sector_speed(test_args: Args) -> None:
    sb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_semiminor_axis(sb, test_args.a, test_args.m)
    with raises(TypeError):
        law.calculate_semiminor_axis(100, test_args.a, test_args.m)


def test_bad_length(test_args: Args) -> None:
    ab = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_semiminor_axis(test_args.s, ab, test_args.m)
    with raises(TypeError):
        law.calculate_semiminor_axis(test_args.s, 100, test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_semiminor_axis(test_args.s, test_args.a, mb)
    with raises(TypeError):
        law.calculate_semiminor_axis(test_args.s, test_args.a, 100)
