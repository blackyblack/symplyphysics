from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal,)

from symplyphysics.laws.electricity.circuits import characteristic_resistance_of_rectangular_waveguide_for_transverse_electric_waves as resistance_law

## The critical wavelength is 17.9 millimeters.The wavelength is 10 millimeters. The characteristic resistance
## of the medium is 254.167 ohms. Then the resistance in the waveguide will be 306.45 ohms.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", ["resistance_of_medium", "wavelength", "critical_wavelength"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    resistance_of_medium = Quantity(254.167 * units.ohm)
    wavelength = Quantity(10 * units.millimeter)
    critical_wavelength = Quantity(17.9 * units.millimeter)
    return Args(resistance_of_medium=resistance_of_medium, wavelength=wavelength, critical_wavelength=critical_wavelength)


def test_basic_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_resistance(test_args.resistance_of_medium, test_args.wavelength,
        test_args.critical_wavelength)
    assert_equal(result, 306.45 * units.ohm)


def test_bad_resistance_of_medium(test_args: Args) -> None:
    bad_resistance_of_medium = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(bad_resistance_of_medium, test_args.wavelength, test_args.critical_wavelength)


def test_bad_wavelength(test_args: Args) -> None:
    bad_wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(test_args.resistance_of_medium, bad_wavelength, test_args.critical_wavelength)
    with raises(TypeError):
        resistance_law.calculate_resistance(test_args.resistance_of_medium, 100, test_args.critical_wavelength)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(test_args.resistance_of_medium, test_args.wavelength, bad_wavelength)
    with raises(TypeError):
        resistance_law.calculate_resistance(test_args.resistance_of_medium, test_args.wavelength, 100)
