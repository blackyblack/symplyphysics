from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import boundary_of_thermalization_zone_of_atomized_atoms_in_magnetron as boundary_law

# Description
## The number of collisions of atomized atom with gas atoms is 6.3. The free path length of atomized atom is 62 millimeter.
## Then boundary of thermalization zone is 39.06 centimeter.

Args = namedtuple("Args", ["number_of_collisions_of_atoms", "free_path_length"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    number_of_collisions_of_atoms = 6.3
    free_path_length = Quantity(62 * units.millimeter)

    return Args(number_of_collisions_of_atoms=number_of_collisions_of_atoms,
        free_path_length=free_path_length)


def test_basic_boundary_of_thermalization_zone(test_args: Args) -> None:
    result = boundary_law.calculate_boundary_of_thermalization_zone(test_args.number_of_collisions_of_atoms, test_args.free_path_length)
    assert_equal(result, 39.06 * units.centimeter)


def test_bad_number_of_collisions_of_atoms(test_args: Args) -> None:
    number_of_collisions_of_atoms = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        boundary_law.calculate_boundary_of_thermalization_zone(number_of_collisions_of_atoms, test_args.free_path_length)


def test_bad_free_path_length(test_args: Args) -> None:
    free_path_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        boundary_law.calculate_boundary_of_thermalization_zone(test_args.number_of_collisions_of_atoms, free_path_length)
    with raises(TypeError):
        boundary_law.calculate_boundary_of_thermalization_zone(test_args.number_of_collisions_of_atoms, 100)
