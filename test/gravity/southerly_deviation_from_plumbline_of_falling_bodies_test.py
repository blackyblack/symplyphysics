from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity import southerly_deviation_from_plumbline_of_falling_bodies as law

Args = namedtuple("Args", "tf tr s lat")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    tf = Quantity(4.5 * units.second)
    tr = Quantity(1 * units.day)
    s = Quantity(1.2 * units.centimeter)
    lat = Quantity(56 * units.degree)
    return Args(tf=tf, tr=tr, s=s, lat=lat)


def test_law(test_args: Args) -> None:
    result = law.calculate_southerly_deviation_from_plumbline(test_args.tf, test_args.tr, test_args.s, test_args.lat)
    assert_equal(result, 1.6 * units.micrometer, tolerance=2e-2)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_southerly_deviation_from_plumbline(tb, test_args.tr, test_args.s, test_args.lat)
    with raises(errors.UnitsError):
        law.calculate_southerly_deviation_from_plumbline(test_args.tf, tb, test_args.s, test_args.lat)
    with raises(TypeError):
        law.calculate_southerly_deviation_from_plumbline(100, test_args.tr, test_args.s, test_args.lat)
    with raises(TypeError):
        law.calculate_southerly_deviation_from_plumbline(test_args.tf, 100, test_args.s, test_args.lat)


def test_bad_length(test_args: Args) -> None:
    sb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_southerly_deviation_from_plumbline(test_args.tf, test_args.tr, sb, test_args.lat)
    with raises(TypeError):
        law.calculate_southerly_deviation_from_plumbline(test_args.tf, test_args.tr, 100, test_args.lat)


def test_bad_angle(test_args: Args) -> None:
    latb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_southerly_deviation_from_plumbline(test_args.tf, test_args.tr, test_args.s, latb)
