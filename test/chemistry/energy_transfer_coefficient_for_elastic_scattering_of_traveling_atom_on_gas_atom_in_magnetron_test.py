from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import energy_transfer_coefficient_for_elastic_scattering_of_traveling_atom_on_gas_atom_in_magnetron as coefficient_law

# Description
## The mass of the traveling titanium atom is 7.95e-26 kilograms. The mass of an argon atom is 6.63e-26 kilogram.
## Then the energy transfer coefficient is 0.496.

Args = namedtuple("Args", ["mass_of_traveling_atom", "mass_of_gas_atom"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mass_of_traveling_atom = Quantity(7.95e-26 * units.kilogram)
    mass_of_gas_atom = Quantity(6.63e-26 * units.kilogram)

    return Args(mass_of_traveling_atom=mass_of_traveling_atom, mass_of_gas_atom=mass_of_gas_atom)


def test_basic_energy_transfer_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_energy_transfer_coefficient(test_args.mass_of_traveling_atom,
        test_args.mass_of_gas_atom)
    assert_equal(result, 0.496)


def test_bad_mass_of_traveling_atom(test_args: Args) -> None:
    mass_of_traveling_atom = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_energy_transfer_coefficient(mass_of_traveling_atom,
            test_args.mass_of_gas_atom)
    with raises(TypeError):
        coefficient_law.calculate_energy_transfer_coefficient(100, test_args.mass_of_gas_atom)


def test_bad_mass_of_gas_atom(test_args: Args) -> None:
    mass_of_gas_atom = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_energy_transfer_coefficient(test_args.mass_of_traveling_atom,
            mass_of_gas_atom)
    with raises(TypeError):
        coefficient_law.calculate_energy_transfer_coefficient(test_args.mass_of_traveling_atom, 100)
