from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import voltage_of_equivalent_diode_for_pentode as voltage_law

# Description
## The distance from the cathode to the anode is 1 centimeter. The distance from the cathode to the first grid is 0.5 centimeter.
## The voltage between the anode and the cathode is 150 volt. The voltage of the first grid is 10 volt. The voltage of the second 
## grid is 20 volt. The voltage of the third grid is 25 volt. The coefficient of direct permeability of the first grid is 2. The
## coefficient of direct permeability of the second grid is 3. The coefficient of direct permeability of the third grid is 1.5.
## Then the voltage between the cathode and the anode of an equivalent diode for a pentode will be equal to 256.64 volt.

Args = namedtuple("Args", ["voltage_of_first_grid", "voltage_of_second_grid", "voltage_of_third_grid", "anode_voltage",
                           "coefficient_direct_permeability_of_first_grid", "coefficient_direct_permeability_of_second_grid",
                           "coefficient_direct_permeability_of_third_grid", "distance_to_anode", "distance_to_first_grid"
                           ])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    voltage_of_first_grid = Quantity(10 * units.volt)
    voltage_of_second_grid = Quantity(20 * units.volt)
    voltage_of_third_grid = Quantity(25 * units.volt)
    anode_voltage = Quantity(150 * units.volt)
    coefficient_direct_permeability_of_first_grid = 2.0
    coefficient_direct_permeability_of_second_grid = 3.0
    coefficient_direct_permeability_of_third_grid = 1.5
    distance_to_anode = Quantity(1 * units.centimeter)
    distance_to_first_grid = Quantity(0.5 * units.centimeter)

    return Args(voltage_of_first_grid=voltage_of_first_grid,
                voltage_of_second_grid=voltage_of_second_grid,
                voltage_of_third_grid=voltage_of_third_grid,
                anode_voltage=anode_voltage,
                coefficient_direct_permeability_of_first_grid=coefficient_direct_permeability_of_first_grid,
                coefficient_direct_permeability_of_second_grid=coefficient_direct_permeability_of_second_grid,
                coefficient_direct_permeability_of_third_grid=coefficient_direct_permeability_of_third_grid,
                distance_to_anode=distance_to_anode,
                distance_to_first_grid=distance_to_first_grid)


def test_basic_voltage_of_equivalent_diode(test_args: Args) -> None:
    result = voltage_law.calculate_voltage_of_equivalent_diode(
        test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)
    ## expected_voltage = (test_args.voltage_of_first_grid + test_args.coefficient_direct_permeability_of_first_grid * test_args.voltage_of_second_grid + test_args.coefficient_direct_permeability_of_first_grid * test_args.coefficient_direct_permeability_of_second_grid * test_args.voltage_of_third_grid + test_args.coefficient_direct_permeability_of_first_grid * test_args.coefficient_direct_permeability_of_second_grid * test_args.coefficient_direct_permeability_of_third_grid * test_args.anode_voltage) / (1 + ((test_args.distance_to_anode / test_args.distance_to_first_grid)**(4 / 3)) * test_args.coefficient_direct_permeability_of_first_grid)
    assert_equal(result, 256.64 * units.volt)


def test_bad_voltage_of_first_grid(test_args: Args) -> None:
    voltage_of_first_grid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(voltage_of_first_grid, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)
    with raises(TypeError):
        voltage_law.calculate_voltage_of_equivalent_diode(100, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)


def test_bad_voltage_of_second_grid(test_args: Args) -> None:
    voltage_of_second_grid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, voltage_of_second_grid, test_args.voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)
    with raises(TypeError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, 100, test_args.voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)


def test_bad_voltage_of_third_grid(test_args: Args) -> None:
    voltage_of_third_grid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)
    with raises(TypeError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, 100, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)


def test_bad_anode_voltage(test_args: Args) -> None:
    anode_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)
    with raises(TypeError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, 100, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)


def test_bad_coefficient_direct_permeability_of_first_grid(test_args: Args) -> None:
    coefficient_direct_permeability_of_first_grid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, test_args.anode_voltage, coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)


def test_bad_coefficient_direct_permeability_of_second_grid(test_args: Args) -> None:
    coefficient_direct_permeability_of_second_grid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)


def test_bad_coefficient_direct_permeability_of_third_grid(test_args: Args) -> None:
    coefficient_direct_permeability_of_third_grid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, test_args.distance_to_first_grid)


def test_bad_distance_to_anode(test_args: Args) -> None:
    distance_to_anode = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, distance_to_anode, test_args.distance_to_first_grid)
    with raises(TypeError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, 100, test_args.distance_to_first_grid)


def test_bad_distance_to_first_grid(test_args: Args) -> None:
    distance_to_first_grid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, distance_to_first_grid)
    with raises(TypeError):
        voltage_law.calculate_voltage_of_equivalent_diode(test_args.voltage_of_first_grid, test_args.voltage_of_second_grid, test_args.voltage_of_third_grid, test_args.anode_voltage, test_args.coefficient_direct_permeability_of_first_grid, test_args.coefficient_direct_permeability_of_second_grid, test_args.coefficient_direct_permeability_of_third_grid, test_args.distance_to_anode, 100)
