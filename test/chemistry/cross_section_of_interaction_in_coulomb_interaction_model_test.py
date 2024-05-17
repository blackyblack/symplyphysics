from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import cross_section_of_interaction_in_coulomb_interaction_model as cross_section_law

# Description
## The ionization energy, expressed in voltage, is 14.69 volt.
## Then the cross-sectional area of the interaction is 2.415e-19 meter^2.

Args = namedtuple("Args", ["ionization_energy"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    ionization_energy = Quantity(14.69 * units.volt)

    return Args(ionization_energy=ionization_energy)


def test_basic_cross_sectional_area_of_interaction(test_args: Args) -> None:
    result = cross_section_law.calculate_cross_sectional_area_of_interaction(
        test_args.ionization_energy)
    assert_equal(result, 2.415e-19 * units.meter**2)


def test_bad_ionization_energy() -> None:
    ionization_energy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(ionization_energy)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(100)
