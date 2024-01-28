from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.optics import linear_magnification_from_distance_to_object_and_distance_to_image as magnification


@fixture(name="test_args")
def test_args_fixture():
    distance_to_image = Quantity(0.175 * units.meter)
    distance_to_object = Quantity(-0.07 * units.meter)
    Args = namedtuple("Args", ["distance_to_image", "distance_to_object"])
    return Args(
        distance_to_image=distance_to_image,
        distance_to_object=distance_to_object,
               )

## The minus sign means that the image is real and inverted.
def test_basic_magnification(test_args):
    result = magnification.calculate_magnification(test_args.distance_to_image, test_args.distance_to_object)
    assert result == approx(-2.5, 0.001)


def test_bad_distance(test_args):
    bad_distance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        magnification.calculate_magnification(bad_distance, test_args.distance_to_object)
    with raises(errors.UnitsError):
        magnification.calculate_magnification(test_args.distance_to_image, bad_distance)
    with raises(TypeError):
        magnification.calculate_magnification(100, test_args.distance_to_object)
    with raises(TypeError):
        magnification.calculate_magnification(test_args.distance_to_image, 100)
    with raises(ValueError):
        magnification.calculate_magnification(test_args.distance_to_image, test_args.distance_to_image)   
