from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors,)
from symplyphysics.laws.astronomy import resolution_of_telescope as resolution_law

# Description
## Let the wavelength be 550 nanometers and the diameter of the lens be 4 millimeters.
## Then the resolution of the telescope will be 1.678-4 radians.
## https://alexandr4784.narod.ru/sdvopdf4/sopgl04_56.pdf

Args = namedtuple("Args", ["wavelength", "lens_diameter"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    wavelength = Quantity(550 * units.nanometer)
    lens_diameter = Quantity(4 * units.millimeter)

    return Args(wavelength=wavelength, lens_diameter=lens_diameter)


def test_basic_resolution(test_args: Args) -> None:
    result = resolution_law.calculate_resolution(test_args.wavelength, test_args.lens_diameter)
    assert_equal(result, 1.678e-4)


def test_bad_wavelength(test_args: Args) -> None:
    wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resolution_law.calculate_resolution(wavelength, test_args.lens_diameter)
    with raises(TypeError):
        resolution_law.calculate_resolution(100, test_args.lens_diameter)


def test_bad_lens_diameter(test_args: Args) -> None:
    lens_diameter = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resolution_law.calculate_resolution(test_args.wavelength, lens_diameter)
    with raises(TypeError):
        resolution_law.calculate_resolution(test_args.wavelength, 100)
