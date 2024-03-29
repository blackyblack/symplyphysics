from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.waveguides import power_carried_by_main_wave_of_rectangular_waveguide as power_law

# Description
## The width and height of the waveguide cross section are 4 centimeter. The wavelength is 10 millimeter.
## The resistance of the material filling the waveguide is 254.167 ohm. The maximum electric field intensity
## is 300 volt per meter. Then the power carried by the waveguide is 140.5 milliwatt.

Args = namedtuple("Args", [
    "waveguide_width", "waveguide_height", "wavelength", "material_resistance",
    "electric_intensity"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    waveguide_width = Quantity(4 * units.centimeter)
    waveguide_height = Quantity(4 * units.centimeter)
    wavelength = Quantity(10 * units.millimeter)
    material_resistance = Quantity(254.167 * units.ohm)
    electric_intensity = Quantity(300 * units.volt / units.meter)

    return Args(waveguide_width=waveguide_width,
        waveguide_height=waveguide_height,
        wavelength=wavelength,
        material_resistance=material_resistance,
        electric_intensity=electric_intensity)


def test_basic_waveguide_power(test_args: Args) -> None:
    result = power_law.calculate_waveguide_power(test_args.waveguide_width,
        test_args.waveguide_height, test_args.wavelength, test_args.material_resistance,
        test_args.electric_intensity)
    assert_equal(result, 140.5 * prefixes.milli * units.watt)


def test_bad_waveguide_width(test_args: Args) -> None:
    waveguide_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_waveguide_power(waveguide_width,
            test_args.waveguide_height, test_args.wavelength, test_args.material_resistance,
            test_args.electric_intensity)


def test_bad_waveguide_height(test_args: Args) -> None:
    waveguide_height = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_waveguide_power(test_args.waveguide_width,
            waveguide_height, test_args.wavelength, test_args.material_resistance,
            test_args.electric_intensity)


def test_bad_wavelength(test_args: Args) -> None:
    wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_waveguide_power(test_args.waveguide_width,
            test_args.waveguide_height, wavelength, test_args.material_resistance,
            test_args.electric_intensity)
    with raises(TypeError):
        power_law.calculate_waveguide_power(test_args.waveguide_width,
            test_args.waveguide_height, 100, test_args.material_resistance,
            test_args.electric_intensity)


def test_bad_material_resistance(test_args: Args) -> None:
    material_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_waveguide_power(test_args.waveguide_width, test_args.waveguide_height, test_args.wavelength, material_resistance, test_args.electric_intensity)
    with raises(TypeError):
        power_law.calculate_waveguide_power(test_args.waveguide_width, test_args.waveguide_height, test_args.wavelength, 100, test_args.electric_intensity)


def test_bad_electric_intensity(test_args: Args) -> None:
    electric_intensity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_waveguide_power(test_args.waveguide_width, test_args.waveguide_height, test_args.wavelength, test_args.material_resistance, electric_intensity)
    with raises(TypeError):
        power_law.calculate_waveguide_power(test_args.waveguide_width, test_args.waveguide_height, test_args.wavelength, test_args.material_resistance, 100)
