from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import number_of_collisions_of_atomized_atom_with_gas_atoms_in_magnetron as number_law

# Description
## The initial energy of the atomized atom is 4 electronvolt. The energy of thermal motion
## in a gas-discharge plasma 0.074 electronvolt. The energy transfer coefficient between the atom and the gas atoms is 0.46.
## Then the number of collisions of atomized atom with gas atoms is 6.475.

Args = namedtuple("Args", ["initial_energy", "energy_of_thermal_motion", "energy_transfer_coefficient"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    initial_energy = Quantity(4 * units.electronvolt)
    energy_of_thermal_motion = Quantity(0.074 * units.electronvolt)
    energy_transfer_coefficient = 0.46

    return Args(initial_energy=initial_energy,
        energy_of_thermal_motion=energy_of_thermal_motion,
        energy_transfer_coefficient=energy_transfer_coefficient)


def test_basic_number_of_collisions_of_atoms(test_args: Args) -> None:
    result = number_law.calculate_number_of_collisions_of_atoms(test_args.initial_energy, test_args.energy_of_thermal_motion, test_args.energy_transfer_coefficient)
    assert_equal(result, 6.475)


def test_bad_initial_energy(test_args: Args) -> None:
    initial_energy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        number_law.calculate_number_of_collisions_of_atoms(initial_energy, test_args.energy_of_thermal_motion, test_args.energy_transfer_coefficient)
    with raises(TypeError):
        number_law.calculate_number_of_collisions_of_atoms(100, test_args.energy_of_thermal_motion, test_args.energy_transfer_coefficient)
    with raises(ValueError):
        number_law.calculate_number_of_collisions_of_atoms(test_args.energy_of_thermal_motion, test_args.initial_energy, test_args.energy_transfer_coefficient)  


def test_bad_energy_of_thermal_motion(test_args: Args) -> None:
    energy_of_thermal_motion = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        number_law.calculate_number_of_collisions_of_atoms(test_args.initial_energy, energy_of_thermal_motion, test_args.energy_transfer_coefficient)
    with raises(TypeError):
        number_law.calculate_number_of_collisions_of_atoms(test_args.initial_energy, 100, test_args.energy_transfer_coefficient)


def test_bad_energy_transfer_coefficient(test_args: Args) -> None:
    energy_transfer_coefficient = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        number_law.calculate_number_of_collisions_of_atoms(test_args.initial_energy, test_args.energy_of_thermal_motion, energy_transfer_coefficient)
    with raises(ValueError):
        number_law.calculate_number_of_collisions_of_atoms(test_args.initial_energy, test_args.energy_of_thermal_motion, 2)
