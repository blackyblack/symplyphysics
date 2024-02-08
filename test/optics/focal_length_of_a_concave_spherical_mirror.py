from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.optics import focal_length_of_a_concave_spherical_mirror as mirror_focus_law

# Test example from https://www.test-uz.ru/solutions_view.php?cat=6&num=15.5


@fixture(name="test_args")
def test_args_fixture():
    r = Quantity(40 * prefixes.centi * units.meters)
    Args = namedtuple("Args", ["r"])
    return Args(r=r)


def test_basic_focus_distance(test_args):
    result = mirror_focus_law.calculate_focus_distance(test_args.r)
    assert_equal(result, 20 * prefixes.centi * units.meters)


def test_bad_wave_length(test_args):
    rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mirror_focus_law.calculate_focus_distance(rb)
    with raises(TypeError):
        mirror_focus_law.calculate_focus_distance(100)
