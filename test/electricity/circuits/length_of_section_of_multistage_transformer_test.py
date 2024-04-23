from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import length_of_section_of_multistage_transformer as length_law

# Description
## The wavelength at the upper frequency is 1 millimeter, the wavelength at the lower frequency is 3 millimeter.
## Then the length of the section will be 0.375 millimeter.

Args = namedtuple("Args", ["wavelength_for_upper_frequency", "wavelength_for_lower_frequency"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    wavelength_for_upper_frequency = Quantity(1 * units.millimeter)
    wavelength_for_lower_frequency = Quantity(3 * units.millimeter)

    return Args(wavelength_for_upper_frequency=wavelength_for_upper_frequency,
        wavelength_for_lower_frequency=wavelength_for_lower_frequency)


def test_basic_length_of_section(test_args: Args) -> None:
    result = length_law.calculate_length_of_section(test_args.wavelength_for_upper_frequency,
        test_args.wavelength_for_lower_frequency)
    assert_equal(result, 0.375 * units.millimeter)


def test_bad_wavelength(test_args: Args) -> None:
    bad_wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        length_law.calculate_length_of_section(bad_wavelength,
            test_args.wavelength_for_lower_frequency)
    with raises(TypeError):
        length_law.calculate_length_of_section(100, test_args.wavelength_for_lower_frequency)
    with raises(errors.UnitsError):
        length_law.calculate_length_of_section(test_args.wavelength_for_upper_frequency,
            bad_wavelength)
    with raises(TypeError):
        length_law.calculate_length_of_section(test_args.wavelength_for_upper_frequency, 100)
