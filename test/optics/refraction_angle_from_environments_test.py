from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
    angle_type,
)
from symplyphysics.laws.optics import refraction_angle_from_environments as refraction_law

# Ray of light transfers from air (refractive coefficient is 1.003) to water (refractive coefficient is 1.333) with incedence angle 30 degrees. What is refraction angle?
## We have online calculator which gives us 22.1 refraction angle. https://matematika-club.ru/kalkulyator-otrazheniya-i-prelomleniya-sveta


@fixture(name="test_args")
def test_args_fixture():
    incedence_angle = Quantity(30 * units.degree)
    incedence_media = 1.003
    refractive_media = 1.333
    Args = namedtuple("Args", ["incedence_media", "incedence_angle", "refractive_media"])
    return Args(incedence_media=incedence_media,
        incedence_angle=incedence_angle,
        refractive_media=refractive_media)


def test_basic_angle(test_args):
    result = refraction_law.calculate_refraction_angle(test_args.incedence_angle,
        test_args.incedence_media, test_args.refractive_media)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, angle_type)
    #HACK: angle quantities are not properly processed by 'convert_to'. Convert their 'scale_factor' instead.
    result_angle = convert_to(result.scale_factor, units.degree).subs(units.degree, 1).evalf(4)
    assert result_angle == approx(22.1, 0.01)


def test_angle_with_number(test_args):
    refraction_law.calculate_refraction_angle(0.5, test_args.incedence_media,
        test_args.refractive_media)
    refraction_law.calculate_refraction_angle(-0.5, test_args.incedence_media,
        test_args.refractive_media)


def test_bad_angle(test_args):
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        refraction_law.calculate_refraction_angle(ab, test_args.incedence_media,
            test_args.refractive_media)
    # Test for large angle
    with raises(AssertionError):
        refraction_law.calculate_refraction_angle(100, test_args.incedence_media,
            test_args.refractive_media)
    with raises(AssertionError):
        refraction_law.calculate_refraction_angle(-100, test_args.incedence_media,
            test_args.refractive_media)
