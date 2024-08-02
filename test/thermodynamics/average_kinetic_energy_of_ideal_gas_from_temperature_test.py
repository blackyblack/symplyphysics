from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import average_kinetic_energy_of_ideal_gas_from_temperature as average_kinetic_energy

Args = namedtuple("Args", ["temperature"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    temperature = Quantity(300 * units.kelvin)
    return Args(temperature=temperature)


def test_basic_average_kinetic_energy(test_args: Args) -> None:
    result = average_kinetic_energy.calculate_average_kinetic_energy(test_args.temperature)
    assert_equal(result, 6.21e-21 * units.joule)


def test_bad_temperature() -> None:
    bad_temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        average_kinetic_energy.calculate_average_kinetic_energy(bad_temperature)
    with raises(TypeError):
        average_kinetic_energy.calculate_average_kinetic_energy(100)
