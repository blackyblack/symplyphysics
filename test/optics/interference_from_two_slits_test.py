from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.optics import \
    interference_from_two_slits as interference_two_holes_law

Args = namedtuple("Args", ["x", "d", "l"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    x = Quantity(5 * prefixes.centi * units.meters)
    d = Quantity(500 * prefixes.micro * units.meters)
    l = Quantity(1.5 * units.meters)
    return Args(x=x, d=d, l=l)


def test_basic_travel_difference(test_args: Args) -> None:
    result = interference_two_holes_law.calculate_travel_difference(test_args.x, test_args.d,
        test_args.l)
    assert_equal(result, 16.67 * prefixes.micro * units.meters)


def test_bad_coordinate(test_args: Args) -> None:
    xb = Quantity(5 * units.coulomb)
    with raises(errors.UnitsError):
        interference_two_holes_law.calculate_travel_difference(xb, test_args.d, test_args.l)
    with raises(TypeError):
        interference_two_holes_law.calculate_travel_difference(100, test_args.d, test_args.l)


def test_bad_distance_to_picture(test_args: Args) -> None:
    lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        interference_two_holes_law.calculate_travel_difference(test_args.x, test_args.d, lb)
    with raises(TypeError):
        interference_two_holes_law.calculate_travel_difference(test_args.x, test_args.d, 100)


def test_bad_distance_between_holes(test_args: Args) -> None:
    db = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        interference_two_holes_law.calculate_travel_difference(test_args.x, db, test_args.l)
    with raises(TypeError):
        interference_two_holes_law.calculate_travel_difference(test_args.x, 100, test_args.l)
