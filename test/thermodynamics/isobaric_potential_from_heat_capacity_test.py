from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.thermodynamics import isobaric_potential_from_heat_capacity as potential_law

# Description
## The decomposition of barium carbonate by the reaction: BaCO3 = BaO + CO2.
## At a temperature of 298 kelvin, the thermal effect of the reaction is 552866 joule per mole, the standart_change_entropy change is 155 [joule / (mole * kelvin)],
## The isobaric heat capacity is -8.79 [joule / (mole * kelvin)].
## Then the isobaric potential is 506676 joule per mole.
## https://studfile.net/preview/5797036/page:4/

Args = namedtuple("Args", ["standart_thermal_effect", "standart_change_entropy", "temperature", "standart_change_heat_capacity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    standart_thermal_effect = Quantity(552866 * units.joule / units.mole)
    standart_change_entropy = Quantity(155 * (units.joule / units.mole / units.kelvin))
    temperature = Quantity(298 * units.kelvin)
    standart_change_heat_capacity = Quantity(-8.79 * (units.joule / units.mole / units.kelvin))

    return Args(standart_thermal_effect=standart_thermal_effect,
        standart_change_entropy=standart_change_entropy,
        temperature=temperature,
        standart_change_heat_capacity=standart_change_heat_capacity)


def test_basic_standart_change_isobaric_potential(test_args: Args) -> None:
    result = potential_law.calculate_standart_change_isobaric_potential(test_args.standart_thermal_effect,
        test_args.standart_change_entropy, test_args.temperature, test_args.standart_change_heat_capacity)
    assert_equal(result, 506676 * units.joule / units.mole)


def test_bad_standart_thermal_effect(test_args: Args) -> None:
    standart_thermal_effect = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_standart_change_isobaric_potential(standart_thermal_effect, test_args.standart_change_entropy,
            test_args.temperature, test_args.standart_change_heat_capacity)
    with raises(TypeError):
        potential_law.calculate_standart_change_isobaric_potential(100, test_args.standart_change_entropy,
            test_args.temperature, test_args.standart_change_heat_capacity)


def test_bad_standart_change_entropy(test_args: Args) -> None:
    standart_change_entropy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_standart_change_isobaric_potential(test_args.standart_thermal_effect, standart_change_entropy,
            test_args.temperature, test_args.standart_change_heat_capacity)
    with raises(TypeError):
        potential_law.calculate_standart_change_isobaric_potential(test_args.standart_thermal_effect, 100,
            test_args.temperature, test_args.standart_change_heat_capacity)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_standart_change_isobaric_potential(test_args.standart_thermal_effect,
            test_args.standart_change_entropy, temperature, test_args.standart_change_heat_capacity)
    with raises(TypeError):
        potential_law.calculate_standart_change_isobaric_potential(test_args.standart_thermal_effect,
            test_args.standart_change_entropy, 100, test_args.standart_change_heat_capacity)


def test_bad_standart_change_heat_capacity(test_args: Args) -> None:
    standart_change_heat_capacity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_standart_change_isobaric_potential(test_args.standart_thermal_effect,
            test_args.standart_change_entropy, test_args.temperature, standart_change_heat_capacity)
    with raises(TypeError):
        potential_law.calculate_standart_change_isobaric_potential(test_args.standart_thermal_effect,
            test_args.standart_change_entropy, test_args.temperature, 100)
