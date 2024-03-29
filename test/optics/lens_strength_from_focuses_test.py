from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.optics import lens_focus_from_object_and_image as lens_law

#We are having thin lens which images object from 0.4m distance to the same 0.4m distance to image. This is only possible if 0.4 is double focus of this lens, so focus should be 0.2m.

Args = namedtuple("Args", ["object_distance", "image_distance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    object_distance = Quantity(0.4 * units.meter)
    image_distance = Quantity(0.4 * units.meter)
    return Args(object_distance=object_distance, image_distance=image_distance)


def test_basic_focus(test_args: Args) -> None:
    result = lens_law.calculate_focus(test_args.object_distance, test_args.image_distance)
    assert_equal(result, 0.2 * units.meter)


def test_bad_distance(test_args: Args) -> None:
    db = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        lens_law.calculate_focus(db, test_args.image_distance)
    with raises(TypeError):
        lens_law.calculate_focus(100, test_args.image_distance)
    with raises(errors.UnitsError):
        lens_law.calculate_focus(test_args.object_distance, db)
    with raises(TypeError):
        lens_law.calculate_focus(test_args.object_distance, 100)
