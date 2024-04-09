from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.waveguides import critical_wavelength_in_rectangular_waveguide as wavelength_law

## Parameters of a rectangular waveguide: width is 2 centimeters, height is 1 cm. The first index is 1, the second index is 1.
## The critical wavelength will be 17.9 millimeters.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", ["first_index", "second_index", "width", "height"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    first_index = 1
    second_index = 1
    width = Quantity(2 * units.centimeter)
    height = Quantity(1 * units.centimeter)
    return Args(first_index=first_index, second_index=second_index, width=width, height=height)


def test_basic_critical_wavelength(test_args: Args) -> None:
    result = wavelength_law.calculate_critical_wavelength(test_args.first_index,
        test_args.second_index, test_args.width, test_args.height)
    assert_equal(result, 17.9 * units.millimeter)


def test_bad_index(test_args: Args) -> None:
    bad_index = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wavelength_law.calculate_critical_wavelength(bad_index, test_args.second_index,
            test_args.width, test_args.height)
    with raises(errors.UnitsError):
        wavelength_law.calculate_critical_wavelength(test_args.first_index, bad_index,
            test_args.width, test_args.height)
    with raises(ValueError):
        wavelength_law.calculate_critical_wavelength(-1, test_args.second_index, test_args.width,
            test_args.height)
    with raises(ValueError):
        wavelength_law.calculate_critical_wavelength(test_args.first_index, -1, test_args.width,
            test_args.height)


def test_bad_width(test_args: Args) -> None:
    width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wavelength_law.calculate_critical_wavelength(test_args.first_index, test_args.second_index,
            width, test_args.height)
    with raises(TypeError):
        wavelength_law.calculate_critical_wavelength(test_args.first_index, test_args.second_index,
            100, test_args.height)


def test_bad_height(test_args: Args) -> None:
    height = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wavelength_law.calculate_critical_wavelength(test_args.first_index, test_args.second_index,
            test_args.width, height)
    with raises(TypeError):
        wavelength_law.calculate_critical_wavelength(test_args.first_index, test_args.second_index,
            test_args.width, 100)
