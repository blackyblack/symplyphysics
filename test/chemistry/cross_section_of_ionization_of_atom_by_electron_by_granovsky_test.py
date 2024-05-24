from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.chemistry import cross_section_of_ionization_of_atom_by_electron_by_granovsky as cross_section_law

## The ionization energy of atom is equal to 15.8 electronvolt. The energy of the ionizing electron is 45 electronvolt.
## The maximum value of the ionization cross-sectional area is 3.1e-20 [meter^2].
## The electron energy corresponding to the maximum value of the ionization cross-sectional area is 110 electronvolt.
## Then the cross-sectional area of the ionization will be 1.916e-20 [meter^2].

Args = namedtuple("Args", [
    "ionization_energy", "energy_of_electron", "maximum_cross_sectional_area_of_ionization", "energy_of_electron_at_max_area",
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    ionization_energy = Quantity(15.8 * units.electronvolt)
    energy_of_electron = Quantity(45 * units.electronvolt)
    maximum_cross_sectional_area_of_ionization = Quantity(3.1e-20 * units.meter**2)
    energy_of_electron_at_max_area = Quantity(110 * units.electronvolt)

    return Args(ionization_energy=ionization_energy,
        energy_of_electron=energy_of_electron,
        maximum_cross_sectional_area_of_ionization=maximum_cross_sectional_area_of_ionization,
        energy_of_electron_at_max_area=energy_of_electron_at_max_area)


def test_basic_cross_sectional_area_of_ionization(test_args: Args) -> None:
    result = cross_section_law.calculate_cross_sectional_area_of_ionization(
        test_args.ionization_energy, test_args.energy_of_electron, test_args.maximum_cross_sectional_area_of_ionization,
        test_args.energy_of_electron_at_max_area)
    assert_equal(result, 1.916E-20 * units.meter**2)


def test_bad_ionization_energy(test_args: Args) -> None:
    ionization_energy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(ionization_energy,
            test_args.energy_of_electron, test_args.maximum_cross_sectional_area_of_ionization,
            test_args.energy_of_electron_at_max_area)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(100,
            test_args.energy_of_electron, test_args.maximum_cross_sectional_area_of_ionization,
            test_args.energy_of_electron_at_max_area)


def test_bad_energy_of_electron(test_args: Args) -> None:
    energy_of_electron = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(test_args.ionization_energy,
            energy_of_electron, test_args.maximum_cross_sectional_area_of_ionization, test_args.energy_of_electron_at_max_area)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(test_args.ionization_energy,
            100, test_args.maximum_cross_sectional_area_of_ionization, test_args.energy_of_electron_at_max_area)


def test_bad_maximum_cross_sectional_area_of_ionization(test_args: Args) -> None:
    maximum_cross_sectional_area_of_ionization = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(test_args.ionization_energy,
            test_args.energy_of_electron, maximum_cross_sectional_area_of_ionization, test_args.energy_of_electron_at_max_area)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(test_args.ionization_energy,
            test_args.energy_of_electron, 100, test_args.energy_of_electron_at_max_area)


def test_bad_energy_of_electron_at_max_area(test_args: Args) -> None:
    energy_of_electron_at_max_area = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(test_args.ionization_energy,
            test_args.energy_of_electron, test_args.maximum_cross_sectional_area_of_ionization, energy_of_electron_at_max_area)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(test_args.ionization_energy,
            test_args.energy_of_electron, test_args.maximum_cross_sectional_area_of_ionization, 100)
