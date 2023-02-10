from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)

from symplyphysics.laws.optics import lens_strength_from_focuses as lens_law

#We are having thin lens which images object from 0.4m distance from object to the same 0.4m distance to image. This is only possible if 0.4 is double focus of this lens, so focus should be 0.2m.

@fixture
def test_args():
    object_distance = units.Quantity('object_distance')
    SI.set_quantity_dimension(object_distance, units.length)
    SI.set_quantity_scale_factor(object_distance, 0.4 * units.meter)
    image_distance = units.Quantity('image_distance')
    SI.set_quantity_dimension(image_distance, units.length)
    SI.set_quantity_scale_factor(image_distance, 0.4 * units.meter)

    Args = namedtuple('Args', ['object_distance', 'image_distance'])
    return Args(object_distance = object_distance, image_distance = image_distance)

def test_basic_focus(test_args):
    result = lens_law.calculate_focus(test_args.object_distance, test_args.image_distance)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_focus = convert_to(result, units.meter).subs(units.meter, 1).evalf(4)
    assert result_focus == approx(0.2, 0.0001)

def test_bad_distance(test_args):
    db = units.Quantity('db')
    SI.set_quantity_dimension(db, units.charge)
    SI.set_quantity_scale_factor(db, 1 * units.coulomb)

    with raises(errors.UnitsError):
        lens_law.calculate_focus(db, test_args.image_distance)

    with raises(TypeError):
        lens_law.calculate_focus(100, test_args.image_distance)

    with raises(errors.UnitsError):
        lens_law.calculate_focus(test_args.object_distance, db)

    with raises(TypeError):
        lens_law.calculate_focus(test_args.object_distance, 100)

