from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to, prefixes,
)
from symplyphysics.laws.optics import optical_distance_from_geometric_distance_and_refractive_index as optical_distance_law

# 1st example: https://mydocx.ru/1-48290.html
# A glass plate with a thickness of h = 1 mm was placed in the path of the beam traveling in the air.
# How much will the optical path length of the beam change if the beam falls on the plate (n = 1.5) normally

# Answer: L = n * h = 1.5 mm

@fixture(name="test_args")
def test_args_fixture():
    h = Quantity(1 * prefixes.milli * units.meters)
    n = 1.5
    Args = namedtuple("Args", ["h", "n"])
    return Args(h=h, n=n)


def test_basic_travel_difference(test_args):
    result = optical_distance_law.calculate_optical_distance(test_args.h, test_args.n)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_travel_difference = convert_to(result, prefixes.milli * units.meters).evalf(3)
    assert result_travel_difference == approx(1.5, 0.1)


def test_bad_distance(test_args):
    hb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        optical_distance_law.calculate_optical_distance(hb, test_args.n)
    with raises(TypeError):
        optical_distance_law.calculate_optical_distance(100, test_args.n)


def test_bad_refractive_index(test_args):
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        optical_distance_law.calculate_optical_distance(test_args.h, nb)
