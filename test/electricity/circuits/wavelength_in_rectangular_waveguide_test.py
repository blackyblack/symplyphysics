from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal,)

from symplyphysics.laws.electricity.circuits import wavelength_in_rectangular_waveguide as wavelength_law

## The critical wavelength is 17.9 millimeters.The wavelength is 10 millimeters.
## Then the wavelength in the waveguide will be 12.057 millimeter.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", ["wavelength", "critical_wavelength"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    wavelength = Quantity(10 * units.millimeter)
    critical_wavelength = Quantity(17.9 * units.millimeter)
    return Args(wavelength=wavelength, critical_wavelength=critical_wavelength)


def test_basic_wavelength_waveguide(test_args: Args) -> None:
    result = wavelength_law.calculate_wavelength_waveguide(test_args.wavelength, test_args.critical_wavelength)
    assert_equal(result, 12.057 * units.millimeter)


def test_bad_wavelength(test_args: Args) -> None:
    bad_wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wavelength_law.calculate_wavelength_waveguide(bad_wavelength, test_args.critical_wavelength)
    with raises(TypeError):
        wavelength_law.calculate_wavelength_waveguide(100, test_args.critical_wavelength)
    with raises(errors.UnitsError):
        wavelength_law.calculate_wavelength_waveguide(test_args.wavelength, bad_wavelength)
    with raises(TypeError):
        wavelength_law.calculate_wavelength_waveguide(test_args.wavelength, 100)
