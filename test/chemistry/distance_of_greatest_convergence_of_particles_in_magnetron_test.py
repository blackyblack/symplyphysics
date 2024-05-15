from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import distance_of_greatest_convergence_of_particles_in_magnetron as distance_law

# Description
## The discharge voltage is 468 volt. The atomic number of the gas atom is 18.
## Then the distance of closest approach is 0.1915 nanometer.

Args = namedtuple("Args", ["discharge_voltage", "number_of_gas_atom"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    discharge_voltage = Quantity(468 * units.volt)
    number_of_gas_atom = 18

    return Args(discharge_voltage=discharge_voltage,
        number_of_gas_atom=number_of_gas_atom)


def test_basic_distance_of_convergence_of_particles(test_args: Args) -> None:
    result = distance_law.calculate_distance_of_convergence_of_particles(test_args.discharge_voltage, test_args.number_of_gas_atom)
    assert_equal(result, 0.1915 * units.nanometer)


def test_bad_discharge_voltage(test_args: Args) -> None:
    discharge_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        distance_law.calculate_distance_of_convergence_of_particles(discharge_voltage, test_args.number_of_gas_atom)
    with raises(TypeError):
        distance_law.calculate_distance_of_convergence_of_particles(100, test_args.number_of_gas_atom)


def test_bad_number_of_gas_atom(test_args: Args) -> None:
    number_of_gas_atom = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        distance_law.calculate_distance_of_convergence_of_particles(test_args.discharge_voltage, number_of_gas_atom)
