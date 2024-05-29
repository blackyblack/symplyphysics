from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.chemistry import concentration_of_ions_in_gas_plasma as concentration_law

# Description
## The probe current is 6.7 kiloampere. The temperature of electrons in plasma is 6000 kelvin.
## The probe surface area is equal to 10 [millimeter^2].
## Then the ion concentration is 3.476e+16 [1 / centimeter^3].

Args = namedtuple("Args", ["probe_current", "electron_temperature", "area_probe_surface"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    probe_current = Quantity(6.7 * prefixes.kilo * units.ampere)
    electron_temperature = Quantity(6000 * units.kelvin)
    area_probe_surface = Quantity(10 * units.millimeter**2)

    return Args(probe_current=probe_current,
        electron_temperature=electron_temperature,
        area_probe_surface=area_probe_surface)


def test_basic_ion_concentration(test_args: Args) -> None:
    result = concentration_law.calculate_ion_concentration(test_args.probe_current, test_args.electron_temperature,
        test_args.area_probe_surface)
    assert_equal(result, 3.476e+16 / units.centimeter**3)


def test_bad_probe_current(test_args: Args) -> None:
    probe_current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        concentration_law.calculate_ion_concentration(probe_current, test_args.electron_temperature,
            test_args.area_probe_surface)
    with raises(TypeError):
        concentration_law.calculate_ion_concentration(100, test_args.electron_temperature,
            test_args.area_probe_surface)


def test_bad_electron_temperature(test_args: Args) -> None:
    electron_temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        concentration_law.calculate_ion_concentration(test_args.probe_current, electron_temperature,
            test_args.area_probe_surface)
    with raises(TypeError):
        concentration_law.calculate_ion_concentration(test_args.probe_current, 100,
            test_args.area_probe_surface)


def test_bad_area_probe_surface(test_args: Args) -> None:
    area_probe_surface = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        concentration_law.calculate_ion_concentration(test_args.probe_current, test_args.electron_temperature,
            area_probe_surface)
    with raises(TypeError):
        concentration_law.calculate_ion_concentration(test_args.probe_current, test_args.electron_temperature, 100)
