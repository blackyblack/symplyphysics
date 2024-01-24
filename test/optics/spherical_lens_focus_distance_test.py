from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to, prefixes,
)
from symplyphysics.laws.optics import spherical_lens_focus_distance as spherical_lens_law

# Test example: https://www.test-uz.ru/solutions_view.php?cat=6&num=15.5


@fixture(name="test_args")
def test_args_fixture():
    r = Quantity(40 * prefixes.centi * units.meters)
    Args = namedtuple("Args", ["r"])
    return Args(r=r)


def test_basic_focus_distance(test_args):
    result = spherical_lens_law.calculate_focus_distance(test_args.r)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_focus = convert_to(result, units.meter).evalf(4)
    assert result_focus == approx(0.2, 0.01)


def test_bad_radius(test_args):
    rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        spherical_lens_law.calculate_focus_distance(rb)
    with raises(TypeError):
        spherical_lens_law.calculate_focus_distance(100)
