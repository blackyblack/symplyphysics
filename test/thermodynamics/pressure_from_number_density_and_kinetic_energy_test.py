from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)

from symplyphysics.laws.thermodynamics import pressure_from_number_density_and_kinetic_energy as ideal_gas_pressure

Args = namedtuple("Args", ["average_kinetic_energy", "molecules_concentration"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    average_kinetic_energy = Quantity(2.4e-20 * units.joule)
    molecules_concentration = Quantity(5e24 * (1 / units.meter**3))
    return Args(average_kinetic_energy=average_kinetic_energy,
        molecules_concentration=molecules_concentration)


def test_basic_pressure(test_args: Args) -> None:
    result = ideal_gas_pressure.calculate_pressure(test_args.molecules_concentration,
        test_args.average_kinetic_energy)
    assert_equal(result, 80000 * units.pascal)


def test_bad_energy(test_args: Args) -> None:
    bad_energy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ideal_gas_pressure.calculate_pressure(test_args.molecules_concentration, bad_energy)
    with raises(TypeError):
        ideal_gas_pressure.calculate_pressure(test_args.molecules_concentration, 100)


def test_bad_concentration(test_args: Args) -> None:
    bad_concentration = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ideal_gas_pressure.calculate_pressure(bad_concentration, test_args.average_kinetic_energy)
    with raises(TypeError):
        ideal_gas_pressure.calculate_pressure(100, test_args.average_kinetic_energy)
