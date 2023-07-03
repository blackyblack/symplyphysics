from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.optics import lens_focus_from_object_and_image as lens_law

#We are having thin lens which images object from 0.4m distance to the same 0.4m distance to image. This is only possible if 0.4 is double focus of this lens, so focus should be 0.2m.


@fixture(name="test_args")
def test_args_fixture():
    object_distance = Quantity(0.4 * units.meter)
    image_distance = Quantity(0.4 * units.meter)
    Args = namedtuple("Args", ["object_distance", "image_distance"])
    return Args(object_distance=object_distance, image_distance=image_distance)


def test_basic_focus(test_args):
    result = lens_law.calculate_focus(test_args.object_distance, test_args.image_distance)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_focus = convert_to(result, units.meter).evalf(4)
    assert result_focus == approx(0.2, 0.0001)


def test_bad_distance(test_args):
    db = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        lens_law.calculate_focus(db, test_args.image_distance)
    with raises(TypeError):
        lens_law.calculate_focus(100, test_args.image_distance)
    with raises(errors.UnitsError):
        lens_law.calculate_focus(test_args.object_distance, db)
    with raises(TypeError):
        lens_law.calculate_focus(test_args.object_distance, 100)
