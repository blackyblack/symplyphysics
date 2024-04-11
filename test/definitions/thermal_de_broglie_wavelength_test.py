from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    quantities,
)
from symplyphysics.definitions import thermal_de_broglie_wavelength

# Description
## The thermal de Broglie wavelength of Argon gas (m = 39.948 amu) at room temperature is 0.16 Å.

Args = namedtuple("Args", "m t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(39.948 * units.amu)
    t = quantities.standard_laboratory_temperature
    return Args(m=m, t=t)


def test_law(test_args: Args) -> None:
    result = thermal_de_broglie_wavelength.calculate_thermal_wavelength(test_args.m, test_args.t)
    assert_equal(result, 0.16 * units.angstrom)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thermal_de_broglie_wavelength.calculate_thermal_wavelength(mb, test_args.t)
    with raises(TypeError):
        thermal_de_broglie_wavelength.calculate_thermal_wavelength(100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thermal_de_broglie_wavelength.calculate_thermal_wavelength(test_args.m, tb)
    with raises(TypeError):
        thermal_de_broglie_wavelength.calculate_thermal_wavelength(test_args.m, 100)
