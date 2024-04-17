from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
    prefixes,
)
from symplyphysics.laws.electricity.waveguides import attenuation_coefficient_in_metal_in_rectangular_waveguide_for_transverse_electric_waves as coefficient_law

# Description
## The surface resistance is 2.576 milliohm.
## The first index and the second index are 1.
## The width is 2 centimeter, the height is 1 centimeter.
## The resistance of medium is 254.167 ohm.
## The critical wavelength is 17.9 millimeters. The wavelength is 10 millimeters.
## Then the attenuation coefficient will be 0.001415 [1 / meter].

Args = namedtuple("Args", [
    "surface_resistance", "first_index", "second_index", "width", "height", "resistance_of_medium",
    "signal_wavelength", "critical_wavelength"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    surface_resistance = Quantity(2.576 * prefixes.milli * units.ohm)
    first_index = 1
    second_index = 1
    width = Quantity(2 * units.centimeter)
    height = Quantity(1 * units.centimeter)

    resistance_of_medium = Quantity(254.167 * units.ohm)
    signal_wavelength = Quantity(10 * units.millimeter)
    critical_wavelength = Quantity(17.9 * units.millimeter)

    return Args(surface_resistance=surface_resistance,
        first_index=first_index,
        second_index=second_index,
        width=width,
        height=height,
        resistance_of_medium=resistance_of_medium,
        signal_wavelength=signal_wavelength,
        critical_wavelength=critical_wavelength)


def test_basic_attenuation_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
        test_args.first_index, test_args.second_index, test_args.width, test_args.height,
        test_args.resistance_of_medium, test_args.signal_wavelength, test_args.critical_wavelength)
    assert_equal(result, 0.001415 * (1 / units.meter))


def test_bad_surface_resistance(test_args: Args) -> None:
    surface_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(surface_resistance, test_args.first_index,
            test_args.second_index, test_args.width, test_args.height,
            test_args.resistance_of_medium, test_args.signal_wavelength,
            test_args.critical_wavelength)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(100, test_args.first_index,
            test_args.second_index, test_args.width, test_args.height,
            test_args.resistance_of_medium, test_args.signal_wavelength,
            test_args.critical_wavelength)


def test_bad_index(test_args: Args) -> None:
    bad_index = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance, bad_index,
            test_args.second_index, test_args.width, test_args.height,
            test_args.resistance_of_medium, test_args.signal_wavelength,
            test_args.critical_wavelength)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, bad_index, test_args.width, test_args.height,
            test_args.resistance_of_medium, test_args.signal_wavelength,
            test_args.critical_wavelength)
    with raises(ValueError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, 0, test_args.width, test_args.height,
            test_args.resistance_of_medium, test_args.signal_wavelength,
            test_args.critical_wavelength)


def test_bad_width(test_args: Args) -> None:
    width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, test_args.second_index, width, test_args.height,
            test_args.resistance_of_medium, test_args.signal_wavelength,
            test_args.critical_wavelength)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, test_args.second_index, 100, test_args.height,
            test_args.resistance_of_medium, test_args.signal_wavelength,
            test_args.critical_wavelength)


def test_bad_height(test_args: Args) -> None:
    height = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, test_args.second_index, test_args.width, height,
            test_args.resistance_of_medium, test_args.signal_wavelength,
            test_args.critical_wavelength)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, test_args.second_index, test_args.width, 100,
            test_args.resistance_of_medium, test_args.signal_wavelength,
            test_args.critical_wavelength)


def test_bad_resistance_of_medium(test_args: Args) -> None:
    resistance_of_medium = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, test_args.second_index, test_args.width, test_args.height,
            resistance_of_medium, test_args.signal_wavelength, test_args.critical_wavelength)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, test_args.second_index, test_args.width, test_args.height, 100,
            test_args.signal_wavelength, test_args.critical_wavelength)


def test_bad_wavelength(test_args: Args) -> None:
    wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, test_args.second_index, test_args.width, test_args.height,
            test_args.resistance_of_medium, wavelength, test_args.critical_wavelength)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, test_args.second_index, test_args.width, test_args.height,
            test_args.resistance_of_medium, 100, test_args.critical_wavelength)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, test_args.second_index, test_args.width, test_args.height,
            test_args.resistance_of_medium, test_args.signal_wavelength, wavelength)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.first_index, test_args.second_index, test_args.width, test_args.height,
            test_args.resistance_of_medium, test_args.signal_wavelength, 100)
