from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.thermodynamics import isobaric_potential_from_heat_capacity as potential_law

# Description
## The decomposition of barium carbonate by the reaction: BaCO3 = BaO + CO2.
## At a temperature of 298 kelvin, the standard thermal effect of the reaction is 552866 joule per mole, the standard change of entropy is 155 [joule / (mole * kelvin)],
## The standard change of isobaric heat capacity is -8.79 [joule / (mole * kelvin)].
## Then the standard change of isobaric potential is 506676 joule per mole.
## https://studfile.net/preview/5797036/page:4/

Args = namedtuple("Args", ["standard_thermal_effect", "standard_change_entropy", "temperature", "standard_change_heat_capacity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    standard_thermal_effect = Quantity(552866 * units.joule / units.mole)
    standard_change_entropy = Quantity(155 * (units.joule / units.mole / units.kelvin))
    temperature = Quantity(298 * units.kelvin)
    standard_change_heat_capacity = Quantity(-8.79 * (units.joule / units.mole / units.kelvin))

    return Args(standard_thermal_effect=standard_thermal_effect,
        standard_change_entropy=standard_change_entropy,
        temperature=temperature,
        standard_change_heat_capacity=standard_change_heat_capacity)


def test_basic_standard_change_isobaric_potential(test_args: Args) -> None:
    result = potential_law.calculate_standard_change_isobaric_potential(test_args.standard_thermal_effect,
        test_args.standard_change_entropy, test_args.temperature, test_args.standard_change_heat_capacity)
    assert_equal(result, 506676 * units.joule / units.mole)


def test_bad_standard_thermal_effect(test_args: Args) -> None:
    standard_thermal_effect = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_standard_change_isobaric_potential(standard_thermal_effect, test_args.standard_change_entropy,
            test_args.temperature, test_args.standard_change_heat_capacity)
    with raises(TypeError):
        potential_law.calculate_standard_change_isobaric_potential(100, test_args.standard_change_entropy,
            test_args.temperature, test_args.standard_change_heat_capacity)


def test_bad_standard_change_entropy(test_args: Args) -> None:
    standard_change_entropy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_standard_change_isobaric_potential(test_args.standard_thermal_effect, standard_change_entropy,
            test_args.temperature, test_args.standard_change_heat_capacity)
    with raises(TypeError):
        potential_law.calculate_standard_change_isobaric_potential(test_args.standard_thermal_effect, 100,
            test_args.temperature, test_args.standard_change_heat_capacity)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_standard_change_isobaric_potential(test_args.standard_thermal_effect,
            test_args.standard_change_entropy, temperature, test_args.standard_change_heat_capacity)
    with raises(TypeError):
        potential_law.calculate_standard_change_isobaric_potential(test_args.standard_thermal_effect,
            test_args.standard_change_entropy, 100, test_args.standard_change_heat_capacity)


def test_bad_standard_change_heat_capacity(test_args: Args) -> None:
    standard_change_heat_capacity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_standard_change_isobaric_potential(test_args.standard_thermal_effect,
            test_args.standard_change_entropy, test_args.temperature, standard_change_heat_capacity)
    with raises(TypeError):
        potential_law.calculate_standard_change_isobaric_potential(test_args.standard_thermal_effect,
            test_args.standard_change_entropy, test_args.temperature, 100)
