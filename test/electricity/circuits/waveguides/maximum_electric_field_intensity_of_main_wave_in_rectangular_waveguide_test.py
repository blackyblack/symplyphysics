from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.circuits.waveguides import maximum_electric_field_intensity_of_main_wave_in_rectangular_waveguide as intensity_law

## The relative permittivity of the material filling the waveguide is 2.2. The width of the cross-section
## of the waveguide is 4 centimeter, the wavelength is 10 millimeter. The magnetic field strength is 100 ampere
## per meter. Then the maximum electric field strength will be 203 kilovolt per meter.

Args = namedtuple("Args",
    ["relative_permittivity", "waveguide_width", "wavelength", "magnetic_intensity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 2.2
    waveguide_width = Quantity(4 * units.centimeter)
    wavelength = Quantity(10 * units.millimeter)
    magnetic_intensity = Quantity(100 * units.ampere / units.meter)
    return Args(relative_permittivity=relative_permittivity,
        waveguide_width=waveguide_width,
        wavelength=wavelength,
        magnetic_intensity=magnetic_intensity)


def test_basic_maximum_electric_intensity(test_args: Args) -> None:
    result = intensity_law.calculate_maximum_electric_intensity(test_args.relative_permittivity,
        test_args.waveguide_width, test_args.wavelength, test_args.magnetic_intensity)
    assert_equal(result, 203 * prefixes.kilo * units.volt / units.meter, relative_tolerance=2e-3)


def test_bad_relative_permittivity(test_args: Args) -> None:
    bad_relative_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_maximum_electric_intensity(bad_relative_permittivity,
            test_args.waveguide_width, test_args.wavelength, test_args.magnetic_intensity)


def test_bad_waveguide_width(test_args: Args) -> None:
    waveguide_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_maximum_electric_intensity(test_args.relative_permittivity,
            waveguide_width, test_args.wavelength, test_args.magnetic_intensity)
    with raises(TypeError):
        intensity_law.calculate_maximum_electric_intensity(test_args.relative_permittivity, 100,
            test_args.wavelength, test_args.magnetic_intensity)


def test_bad_wavelength(test_args: Args) -> None:
    wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_maximum_electric_intensity(test_args.relative_permittivity,
            test_args.waveguide_width, wavelength, test_args.magnetic_intensity)
    with raises(TypeError):
        intensity_law.calculate_maximum_electric_intensity(test_args.relative_permittivity,
            test_args.waveguide_width, 100, test_args.magnetic_intensity)


def test_bad_magnetic_intensity(test_args: Args) -> None:
    magnetic_intensity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_maximum_electric_intensity(test_args.relative_permittivity,
            test_args.waveguide_width, test_args.wavelength, magnetic_intensity)
    with raises(TypeError):
        intensity_law.calculate_maximum_electric_intensity(test_args.relative_permittivity,
            test_args.waveguide_width, test_args.wavelength, 100)
