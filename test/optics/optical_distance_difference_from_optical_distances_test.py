from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.optics import optical_distance_difference_from_optical_distances as optical_difference_law

# 1st example: https://mydocx.ru/1-48290.html
# A glass plate with a thickness of h = 1 mm was placed in the path of the beam traveling in the air.
# How much will the optical path length of the beam change if the beam falls on the plate (n = 1.5) normally

# Answer: delta = L2 - L1 = n * h - 1 * h= 0.5 mm

Args = namedtuple("Args", ["ol1", "ol2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    ol1 = Quantity(1 * prefixes.milli * units.meters)
    ol2 = Quantity(1.5 * prefixes.milli * units.meters)
    return Args(ol1=ol1, ol2=ol2)


def test_basic_distance_difference(test_args: Args) -> None:
    result = optical_difference_law.calculate_optical_difference_distance(
        test_args.ol1, test_args.ol2)
    assert_equal(result, 0.5 * prefixes.milli * units.meters)


def test_bad_distances(test_args: Args) -> None:
    olb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        optical_difference_law.calculate_optical_difference_distance(olb, test_args.ol2)
    with raises(TypeError):
        optical_difference_law.calculate_optical_difference_distance(100, test_args.ol2)

    with raises(errors.UnitsError):
        optical_difference_law.calculate_optical_difference_distance(test_args.ol1, olb)
    with raises(TypeError):
        optical_difference_law.calculate_optical_difference_distance(test_args.ol1, 100)
