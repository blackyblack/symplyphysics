from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.optics import linear_magnification_from_object_height_and_image_height as magnification


@fixture(name="test_args")
def test_args_fixture():
    image_height = Quantity(5 * units.meter)
    object_height = Quantity(3.5 * units.meter)
    Args = namedtuple("Args", ["image_height", "object_height"])
    return Args(
        image_height=image_height,
        object_height=object_height,
    )


def test_basic_magnification(test_args):
    result = magnification.calculate_magnification(test_args.image_height, test_args.object_height)
    assert_approx(result, 1.428)


def test_bad_height(test_args):
    bad_height = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        magnification.calculate_magnification(bad_height, test_args.object_height)
    with raises(errors.UnitsError):
        magnification.calculate_magnification(test_args.image_height, bad_height)
    with raises(TypeError):
        magnification.calculate_magnification(100, test_args.object_height)
    with raises(TypeError):
        magnification.calculate_magnification(test_args.image_height, 100)
