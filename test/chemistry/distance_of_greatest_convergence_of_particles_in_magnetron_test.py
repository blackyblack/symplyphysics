from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import distance_of_greatest_convergence_of_particles_in_magnetron as distance_law

# Description
## The discharge voltage is 468 volt. The atomic number of the first atom is 22.
## The atomic number of the second atom is 18. Then the distance of closest approach is 0.2 nanometer.

Args = namedtuple("Args",
    ["discharge_voltage", "atomic_number_of_first_atom", "atomic_number_of_second_atom"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    discharge_voltage = Quantity(468 * units.volt)
    atomic_number_of_first_atom = 22
    atomic_number_of_second_atom = 18

    return Args(discharge_voltage=discharge_voltage,
        atomic_number_of_first_atom=atomic_number_of_first_atom,
        atomic_number_of_second_atom=atomic_number_of_second_atom)


def test_basic_distance_of_convergence_of_particles(test_args: Args) -> None:
    result = distance_law.calculate_distance_of_convergence_of_particles(
        test_args.discharge_voltage, test_args.atomic_number_of_first_atom,
        test_args.atomic_number_of_second_atom)
    assert_equal(result, 0.2 * units.nanometer, relative_tolerance=0.01)


def test_bad_discharge_voltage(test_args: Args) -> None:
    discharge_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        distance_law.calculate_distance_of_convergence_of_particles(
            discharge_voltage, test_args.atomic_number_of_first_atom,
            test_args.atomic_number_of_second_atom)
    with raises(TypeError):
        distance_law.calculate_distance_of_convergence_of_particles(
            100, test_args.atomic_number_of_first_atom, test_args.atomic_number_of_second_atom)


def test_bad_atomic_numbers(test_args: Args) -> None:
    bad_atomic_number = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        distance_law.calculate_distance_of_convergence_of_particles(
            test_args.discharge_voltage, bad_atomic_number, test_args.atomic_number_of_second_atom)
    with raises(errors.UnitsError):
        distance_law.calculate_distance_of_convergence_of_particles(
            test_args.discharge_voltage, test_args.atomic_number_of_first_atom, bad_atomic_number)
