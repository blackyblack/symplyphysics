from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)

from symplyphysics.core import base_law
from symplyphysics.laws.optics import linear_magnification_from_object_height_and_image_height as magnification_law

Args = namedtuple("Args", ["ho", "hi"])
law = base_law.BaseLaw(magnification_law.law)

@fixture(name="test_args")
def test_args_fixture() -> Args:
    ho = Quantity(7 * units.meters)
    hi = Quantity(5 * units.meters)
    return Args(ho=ho, hi=hi)


def test_basic_linear_magnification(test_args: Args) -> None:
    dict_for_calculating = {
        magnification_law.object_height: test_args.ho,
        magnification_law.image_height: test_args.hi
    }
    result = law.calculate_symbol_value(magnification_law.magnification, dict_for_calculating)
    assert_equal(result, 0.714)


def test_bad_height_object(test_args: Args) -> None:
    hob = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        dict_for_calculating = {
            magnification_law.object_height: hob,
            magnification_law.image_height: test_args.hi
        }
        law.calculate_symbol_value(magnification_law.magnification, dict_for_calculating)
    with raises(TypeError):
        dict_for_calculating = {
            magnification_law.object_height: 100,
            magnification_law.image_height: test_args.hi
        }
        law.calculate_symbol_value(magnification_law.magnification, dict_for_calculating)


def test_bad_height_image(test_args: Args) -> None:
    hib = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        dict_for_calculating = {
            magnification_law.object_height: test_args.ho,
            magnification_law.image_height: hib
        }
        law.calculate_symbol_value(magnification_law.magnification, dict_for_calculating)
    with raises(TypeError):
        dict_for_calculating = {
            magnification_law.object_height: test_args.ho,
            magnification_law.image_height: 100
        }
        law.calculate_symbol_value(magnification_law.magnification, dict_for_calculating)
