from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import space_charge_density_in_diode_device as density_law

# Description
## The voltage on the grid is 100 volt. The distance between the grid and the cathode is 2 centimeter.
## Then the space charge density in the diode device will be 9.84e-7 [coulomb / meter^3].

Args = namedtuple("Args", ["voltage_on_grid", "distance_to_cathode"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    voltage_on_grid = Quantity(100 * units.volt)
    distance_to_cathode = Quantity(2 * units.centimeter)

    return Args(voltage_on_grid=voltage_on_grid, distance_to_cathode=distance_to_cathode)


def test_basic_space_charge_density(test_args: Args) -> None:
    result = density_law.calculate_space_charge_density(test_args.voltage_on_grid,
        test_args.distance_to_cathode)
    assert_equal(result, 9.84e-7 * units.coulomb / units.meter**3)


def test_bad_voltage_on_grid(test_args: Args) -> None:
    voltage_on_grid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        density_law.calculate_space_charge_density(voltage_on_grid, test_args.distance_to_cathode)
    with raises(TypeError):
        density_law.calculate_space_charge_density(100, test_args.distance_to_cathode)


def test_bad_distance_to_cathode(test_args: Args) -> None:
    distance_to_cathode = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        density_law.calculate_space_charge_density(test_args.voltage_on_grid, distance_to_cathode)
    with raises(TypeError):
        density_law.calculate_space_charge_density(test_args.voltage_on_grid, 100)
