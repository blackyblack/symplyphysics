from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import voltage_of_equivalent_diode_for_triode as voltage_law

# Description
## The distance from the cathode to the anode is 4 centimeter. The distance from the cathode to the grid is 1.5 centimeter.
## The voltage between the anode and the cathode is 17 volt. The voltage triode gain is 3. The grid voltage is -2 volt.
## Then the voltage between the cathode and the anode of an equivalent diode for a triode will be equal to 1.642 volt.

Args = namedtuple("Args", [
    "distance_to_anode", "distance_to_grid", "anode_voltage", "voltage_triode_gain", "grid_voltage"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    distance_to_anode = Quantity(4 * units.centimeter)
    distance_to_grid = Quantity(1.5 * units.centimeter)
    anode_voltage = Quantity(17 * units.volt)
    voltage_triode_gain = 3
    grid_voltage = Quantity(-2 * units.volt)

    return Args(distance_to_anode=distance_to_anode,
        distance_to_grid=distance_to_grid,
        anode_voltage=anode_voltage,
        voltage_triode_gain=voltage_triode_gain,
        grid_voltage=grid_voltage)


def test_basic_voltage_of_equivalent_diode(test_args: Args) -> None:
    result = voltage_law.calculate_voltage_of_equivalent_diode(test_args.distance_to_anode,
        test_args.distance_to_grid, test_args.anode_voltage, test_args.voltage_triode_gain,
        test_args.grid_voltage)
    assert_equal(result, 1.642 * units.volt)


def test_bad_distance_to_anode(test_args: Args) -> None:
    distance_to_anode = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(distance_to_anode,
            test_args.distance_to_grid, test_args.anode_voltage, test_args.voltage_triode_gain,
            test_args.grid_voltage)
    with raises(TypeError):
        voltage_law.calculate_voltage_of_equivalent_diode(100, test_args.distance_to_grid,
            test_args.anode_voltage, test_args.voltage_triode_gain, test_args.grid_voltage)


def test_bad_distance_to_grid(test_args: Args) -> None:
    distance_to_grid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.distance_to_anode,
            distance_to_grid, test_args.anode_voltage, test_args.voltage_triode_gain,
            test_args.grid_voltage)
    with raises(TypeError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.distance_to_anode, 100,
            test_args.anode_voltage, test_args.voltage_triode_gain, test_args.grid_voltage)


def test_bad_anode_voltage(test_args: Args) -> None:
    anode_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.distance_to_anode,
            test_args.distance_to_grid, anode_voltage, test_args.voltage_triode_gain,
            test_args.grid_voltage)
    with raises(TypeError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.distance_to_anode,
            test_args.distance_to_grid, 100, test_args.voltage_triode_gain, test_args.grid_voltage)


def test_bad_voltage_triode_gain(test_args: Args) -> None:
    voltage_triode_gain = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.distance_to_anode,
            test_args.distance_to_grid, test_args.anode_voltage, voltage_triode_gain,
            test_args.grid_voltage)


def test_bad_grid_voltage(test_args: Args) -> None:
    grid_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.distance_to_anode,
            test_args.distance_to_grid, test_args.anode_voltage, test_args.voltage_triode_gain,
            grid_voltage)
    with raises(TypeError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.distance_to_anode,
            test_args.distance_to_grid, test_args.anode_voltage, test_args.voltage_triode_gain, 100)
