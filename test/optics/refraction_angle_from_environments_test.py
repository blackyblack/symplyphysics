from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.optics import refraction_angle_from_environments as refraction_law

# Ray of light transfers from air (refractive coefficient is 1.003) to water (refractive coefficient is 1.333) with incidence angle 30 degrees. What is refraction angle?
## We have online calculator which gives us 22.1 refraction angle. https://matematika-club.ru/kalkulyator-otrazheniya-i-prelomleniya-sveta


@fixture(name="test_args")
def test_args_fixture():
    incidence_angle = Quantity(30 * units.degree)
    incidence_media = 1.003
    refractive_media = 1.333
    Args = namedtuple("Args", ["incidence_media", "incidence_angle", "refractive_media"])
    return Args(incidence_media=incidence_media,
        incidence_angle=incidence_angle,
        refractive_media=refractive_media)


def test_basic_angle(test_args):
    result = refraction_law.calculate_refraction_angle(test_args.incidence_angle,
        test_args.incidence_media, test_args.refractive_media)
    assert_equal(result, 22.1 * units.degree)


def test_angle_with_number(test_args):
    refraction_law.calculate_refraction_angle(0.5, test_args.incidence_media,
        test_args.refractive_media)
    refraction_law.calculate_refraction_angle(-0.5, test_args.incidence_media,
        test_args.refractive_media)


def test_bad_angle(test_args):
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        refraction_law.calculate_refraction_angle(ab, test_args.incidence_media,
            test_args.refractive_media)
    # Test for large angle
    with raises(AssertionError):
        refraction_law.calculate_refraction_angle(100, test_args.incidence_media,
            test_args.refractive_media)
    with raises(AssertionError):
        refraction_law.calculate_refraction_angle(-100, test_args.incidence_media,
            test_args.refractive_media)
