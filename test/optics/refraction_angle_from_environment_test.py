from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)

from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

from symplyphysics.laws.optics import refraction_angle_from_environments as refraction_law

# Ray of light transfers from air (reflective coefficient is 1.003) to water (reflective coefficient is 1.333) with incedence angle 30 degrees. What is refraction angle?
## We have online calculator which gives us 22.1 refraction angle. https://matematika-club.ru/kalkulyator-otrazheniya-i-prelomleniya-sveta

@fixture
def test_args():
    incedence_angle = units.Quantity('incedence_angle')
    SI.set_quantity_dimension(incedence_angle, angle_type)
    SI.set_quantity_scale_factor(incedence_angle, 30 * units.degree)

    incedence_media = 1.003
    refractive_media = 1.333
    
    Args = namedtuple('Args', ['incedence_media', 'incedence_angle', 'refractive_media'])
    return Args(incedence_media = incedence_media, incedence_angle = incedence_angle, refractive_media = refractive_media)

def test_basic_angle(test_args):
    result = refraction_law.calculate_angle(test_args.incedence_media, test_args.incedence_angle, test_args.refractive_media)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, angle_type)
    result_angle = convert_to(result, units.degree).subs(units.degree, 1).evalf(4)
    assert result_angle == approx(22.1, 0.01)

def test_bad_angle(test_args):
    ab = units.Quantity('ab')
    SI.set_quantity_dimension(ab, units.charge)
    SI.set_quantity_scale_factor(ab, 1 * units.coulomb)

    with raises(errors.UnitsError):
        refraction_law.calculate_angle(test_args.incedence_media, ab, test_args.refractive_media)

    with raises(TypeError):
        refraction_law.calculate_angle(test_args.incedence_media, 100, test_args.refractive_media)
