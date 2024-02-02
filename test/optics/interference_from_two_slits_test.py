from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to, prefixes,
)
from symplyphysics.laws.optics import \
    interference_from_two_slits as interference_two_holes_law


@fixture(name="test_args")
def test_args_fixture():
    x = Quantity(5 * prefixes.centi * units.meters)
    d = Quantity(500 * prefixes.micro * units.meters)
    l = Quantity(1.5 * units.meters)
    Args = namedtuple("Args", ["x", "d", "l"])
    return Args(x=x, d=d, l=l)


def test_basic_travel_difference(test_args):
    result = interference_two_holes_law.calculate_travel_difference(test_args.x,
                                                                    test_args.d,
                                                                    test_args.l)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_travel_difference = convert_to(result, prefixes.micro * units.meters).evalf(5)
    assert result_travel_difference == approx(16.67, 0.001)


def test_bad_coordinate(test_args):
    xb = Quantity(5 * units.coulomb)
    with raises(errors.UnitsError):
        interference_two_holes_law.calculate_travel_difference(xb, test_args.d,
                                                               test_args.l)
    with raises(TypeError):
        interference_two_holes_law.calculate_travel_difference(100, test_args.d,
                                                               test_args.l)


def test_bad_distance_to_picture(test_args):
    lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        interference_two_holes_law.calculate_travel_difference(test_args.x,
                                                               test_args.d, lb)
    with raises(TypeError):
        interference_two_holes_law.calculate_travel_difference(test_args.x,
                                                               test_args.d, 100)


def test_bad_distance_between_holes(test_args):
    db = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        interference_two_holes_law.calculate_travel_difference(test_args.x, db,
                                                               test_args.l)
    with raises(TypeError):
        interference_two_holes_law.calculate_travel_difference(test_args.x, 100,
                                                               test_args.l)
