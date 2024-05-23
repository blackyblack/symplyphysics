from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.chemistry import cross_section_of_ionization_of_atom_by_electron_by_lotz_drevin as cross_section_law

## The ionization energy of atom is equal to 15.8 electronvolt. The energy of the ionizing electron is 110 electronvolt.
## The first and second calculation coefficients are 1.9 and 0.4. The number of equivalent electrons on the outer shell of an atom is 6.
## Then the cross-sectional area of the ionization will be 3.03e-20 [meter^2].

Args = namedtuple("Args", [
    "ionization_energy", "energy_of_electron", "first_calculation_coefficient", "second_calculation_coefficient", "number_of_electrons"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    ionization_energy = Quantity(15.8 * units.electronvolt)
    energy_of_electron = Quantity(110 * units.electronvolt)
    first_calculation_coefficient = 1.9
    second_calculation_coefficient = 0.4
    number_of_electrons = 6

    return Args(ionization_energy=ionization_energy,
        energy_of_electron=energy_of_electron,
        first_calculation_coefficient=first_calculation_coefficient,
        second_calculation_coefficient=second_calculation_coefficient,
        number_of_electrons=number_of_electrons)


def test_basic_cross_sectional_area_of_ionization(test_args: Args) -> None:
    result = cross_section_law.calculate_cross_sectional_area_of_ionization(
        test_args.ionization_energy, test_args.energy_of_electron, test_args.first_calculation_coefficient,
        test_args.second_calculation_coefficient, test_args.number_of_electrons,)
    assert_equal(result, 3.03e-20 * units.meter**2)


def test_bad_ionization_energy(test_args: Args) -> None:
    ionization_energy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(ionization_energy,
            test_args.energy_of_electron, test_args.first_calculation_coefficient,
            test_args.second_calculation_coefficient, test_args.number_of_electrons)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(100,
            test_args.energy_of_electron, test_args.first_calculation_coefficient,
            test_args.second_calculation_coefficient, test_args.number_of_electrons)


def test_bad_energy_of_electron(test_args: Args) -> None:
    energy_of_electron = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(test_args.ionization_energy,
            energy_of_electron, test_args.first_calculation_coefficient, test_args.second_calculation_coefficient, test_args.number_of_electrons)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(test_args.ionization_energy,
            100, test_args.first_calculation_coefficient, test_args.second_calculation_coefficient, test_args.number_of_electrons)


def test_bad_calculation_coefficients(test_args: Args) -> None:
    calculation_coefficient = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(test_args.ionization_energy,
            test_args.energy_of_electron, calculation_coefficient, test_args.second_calculation_coefficient, test_args.number_of_electrons)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(test_args.ionization_energy,
            test_args.energy_of_electron, test_args.first_calculation_coefficient, calculation_coefficient, test_args.number_of_electrons)


def test_bad_number_of_electrons(test_args: Args) -> None:
    number_of_electrons = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_ionization(test_args.ionization_energy,
            test_args.energy_of_electron, test_args.first_calculation_coefficient, test_args.second_calculation_coefficient, number_of_electrons)
