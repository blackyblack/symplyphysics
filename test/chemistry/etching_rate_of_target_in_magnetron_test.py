from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import etching_rate_of_target_in_magnetron as rate_law

# Description
## The current density of ions incident on the target is 97.566 [ampere / meter^2]. The molar mass of the target atom
## is 47.867 gram per mole. The atomization coefficient of the target atoms is 1.873. The target density is
## 6020 [kilogram / meter^3]. Then the etching rate of the target is 15.06 nanometer per second.

Args = namedtuple("Args", [
    "ion_current_density", "molar_mass_of_target_atom", "sputtering_coefficient", "target_density"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    ion_current_density = Quantity(97.566 * units.ampere / units.meter**2)
    molar_mass_of_target_atom = Quantity(47.867 * units.gram / units.mol)
    sputtering_coefficient = 1.873
    target_density = Quantity(6020 * units.kilogram / units.meter**3)

    return Args(ion_current_density=ion_current_density,
        molar_mass_of_target_atom=molar_mass_of_target_atom,
        sputtering_coefficient=sputtering_coefficient,
        target_density=target_density)


def test_basic_etching_rate(test_args: Args) -> None:
    result = rate_law.calculate_etching_rate(test_args.ion_current_density,
        test_args.molar_mass_of_target_atom, test_args.sputtering_coefficient,
        test_args.target_density)
    assert_equal(result, 15.06 * (units.nanometer / units.second))


def test_bad_ion_current_density(test_args: Args) -> None:
    ion_current_density = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        rate_law.calculate_etching_rate(ion_current_density, test_args.molar_mass_of_target_atom,
            test_args.sputtering_coefficient, test_args.target_density)
    with raises(TypeError):
        rate_law.calculate_etching_rate(100, test_args.molar_mass_of_target_atom,
            test_args.sputtering_coefficient, test_args.target_density)


def test_bad_molar_mass_of_target_atom(test_args: Args) -> None:
    molar_mass_of_target_atom = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        rate_law.calculate_etching_rate(test_args.ion_current_density, molar_mass_of_target_atom,
            test_args.sputtering_coefficient, test_args.target_density)
    with raises(TypeError):
        rate_law.calculate_etching_rate(test_args.ion_current_density, 100,
            test_args.sputtering_coefficient, test_args.target_density)


def test_bad_sputtering_coefficient(test_args: Args) -> None:
    sputtering_coefficient = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        rate_law.calculate_etching_rate(test_args.ion_current_density,
            test_args.molar_mass_of_target_atom, sputtering_coefficient, test_args.target_density)


def test_bad_target_density(test_args: Args) -> None:
    target_density = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        rate_law.calculate_etching_rate(test_args.ion_current_density,
            test_args.molar_mass_of_target_atom, test_args.sputtering_coefficient, target_density)
    with raises(TypeError):
        rate_law.calculate_etching_rate(
            test_args.ion_current_density,
            test_args.molar_mass_of_target_atom,
            test_args.sputtering_coefficient,
            100,
        )
