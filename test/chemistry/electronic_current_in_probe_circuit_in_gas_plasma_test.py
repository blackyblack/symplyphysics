from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.chemistry import electronic_current_in_probe_circuit_in_gas_plasma as current_law

## The probe surface area is equal to 10 [millimeter^2]. The concentration of electrons in plasma is 3.5e16 [1 / centimeter^3].
## The temperature of electrons in plasma is 6000 kelvin. The floating plasma potential is 100 volt. The potential at the location of the probe is 95 volt.
## Then the probe current will be 425.84 milliampere.

Args = namedtuple("Args", [
    "area_probe_surface", "electron_concentration", "plasma_temperature",
    "floating_plasma_potential", "probe_potential"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    area_probe_surface = Quantity(10 * units.millimeter**2)
    electron_concentration = Quantity(3.5e16 / units.centimeter**3)
    plasma_temperature = Quantity(6000 * units.kelvin)
    floating_plasma_potential = Quantity(100 * units.volt)
    probe_potential = Quantity(95 * units.volt)

    return Args(area_probe_surface=area_probe_surface,
        electron_concentration=electron_concentration,
        plasma_temperature=plasma_temperature,
        floating_plasma_potential=floating_plasma_potential,
        probe_potential=probe_potential)


def test_basic_current(test_args: Args) -> None:
    result = current_law.calculate_current(
        test_args.area_probe_surface,
        test_args.electron_concentration,
        test_args.plasma_temperature,
        test_args.floating_plasma_potential,
        test_args.probe_potential,
    )
    assert_equal(result, 425.84 * prefixes.milli * units.ampere)


def test_bad_area_probe_surface(test_args: Args) -> None:
    area_probe_surface = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(area_probe_surface, test_args.electron_concentration,
            test_args.plasma_temperature, test_args.floating_plasma_potential,
            test_args.probe_potential)
    with raises(TypeError):
        current_law.calculate_current(100, test_args.electron_concentration,
            test_args.plasma_temperature, test_args.floating_plasma_potential,
            test_args.probe_potential)


def test_bad_electron_concentration(test_args: Args) -> None:
    electron_concentration = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(test_args.area_probe_surface, electron_concentration,
            test_args.plasma_temperature, test_args.floating_plasma_potential,
            test_args.probe_potential)
    with raises(TypeError):
        current_law.calculate_current(test_args.area_probe_surface, 100,
            test_args.plasma_temperature, test_args.floating_plasma_potential,
            test_args.probe_potential)


def test_bad_plasma_temperature(test_args: Args) -> None:
    plasma_temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(test_args.area_probe_surface,
            test_args.electron_concentration, plasma_temperature,
            test_args.floating_plasma_potential, test_args.probe_potential)
    with raises(TypeError):
        current_law.calculate_current(test_args.area_probe_surface,
            test_args.electron_concentration, 100, test_args.floating_plasma_potential,
            test_args.probe_potential)


def test_bad_floating_plasma_potential(test_args: Args) -> None:
    floating_plasma_potential = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(test_args.area_probe_surface,
            test_args.electron_concentration, test_args.plasma_temperature,
            floating_plasma_potential, test_args.probe_potential)
    floating_plasma_potential = Quantity(1 * units.coulomb)
    with raises(TypeError):
        current_law.calculate_current(test_args.area_probe_surface,
            test_args.electron_concentration, test_args.plasma_temperature, 100,
            test_args.probe_potential)


def test_bad_probe_potential(test_args: Args) -> None:
    probe_potential = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(test_args.area_probe_surface,
            test_args.electron_concentration, test_args.plasma_temperature,
            test_args.floating_plasma_potential, probe_potential)
    with raises(TypeError):
        current_law.calculate_current(test_args.area_probe_surface,
            test_args.electron_concentration, test_args.plasma_temperature,
            test_args.floating_plasma_potential, 100)
