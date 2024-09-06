from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import energy_density_via_permittivity_and_electric_field as energy_density_law

# Description
## It is known that with a permittivity equal to 5 and an electric field intensity equal to 5 Volt / meter,
## energy density of the electric field is (1.106e-10 * 5) Joule / meter**3.
## https://byjus.com/energy-density-formula/
## Example above calculates energy density with vacuum permittivity. We should multiply the result to relative
## permittivity to obtain the medium permittivity.

Args = namedtuple("Args", ["absolute_permittivity", "electric_intensity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    _relative_permittivity = 5
    absolute_permittivity = Quantity(units.vacuum_permittivity * _relative_permittivity)
    electric_intensity = Quantity(5 * (units.volt / units.meter))
    return Args(absolute_permittivity=absolute_permittivity, electric_intensity=electric_intensity)


def test_basic_energy_density(test_args: Args) -> None:
    result = energy_density_law.calculate_energy_density(test_args.absolute_permittivity,
        test_args.electric_intensity)
    assert_equal(result, 1.106e-10 * 5 * units.joule / units.meter**3)


def test_bad_absolute_permittivity(test_args: Args) -> None:
    absolute_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_density_law.calculate_energy_density(absolute_permittivity,
            test_args.electric_intensity)
    with raises(TypeError):
        energy_density_law.calculate_energy_density(100, test_args.electric_intensity)


def test_bad_electric_intensity(test_args: Args) -> None:
    electric_intensity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_density_law.calculate_energy_density(test_args.absolute_permittivity,
            electric_intensity)
    with raises(TypeError):
        energy_density_law.calculate_energy_density(test_args.absolute_permittivity, 100)
