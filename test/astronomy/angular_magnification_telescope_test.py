from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import angular_magnification_telescope as magnification_law

# Description
## With a focal length of the lens equal to 20 meters and a focal length of the eyepiece equal to 0.5 meters,
## the angular magnification of the telescope will be 40.

Args = namedtuple("Args", ["focal_length_lens", "focal_length_eyepiece"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    focal_length_lens = Quantity(20 * units.meter)
    focal_length_eyepiece = Quantity(0.5 * units.meter)

    return Args(focal_length_lens=focal_length_lens, focal_length_eyepiece=focal_length_eyepiece)


def test_basic_angular_magnification(test_args: Args) -> None:
    result = magnification_law.calculate_angular_magnification(test_args.focal_length_lens,
        test_args.focal_length_eyepiece)
    assert_equal(result, 40)


def test_bad_focal_lengths(test_args: Args) -> None:
    bad_focal_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        magnification_law.calculate_angular_magnification(bad_focal_length,
            test_args.focal_length_eyepiece)
    with raises(TypeError):
        magnification_law.calculate_angular_magnification(100, test_args.focal_length_eyepiece)
    with raises(errors.UnitsError):
        magnification_law.calculate_angular_magnification(test_args.focal_length_lens,
            bad_focal_length)
    with raises(TypeError):
        magnification_law.calculate_angular_magnification(test_args.focal_length_lens, 100)
