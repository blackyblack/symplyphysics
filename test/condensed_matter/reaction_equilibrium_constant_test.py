from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, dimensionless)
from symplyphysics.laws.condensed_matter import reaction_equilibrium_constant as constant_law

# Description
## The decomposition of barium carbonate by the reaction: BaCO3 = BaO + CO2.
## At a temperature of 2000 kelvin, the change of isobaric potential is -74150 joule per mole.
## Then equilibrium constant of reaction is 86.409.
## https://studfile.net/preview/5797036/page:4/

Args = namedtuple("Args", ["standart_change_isobaric_potential", "temperature"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    standart_change_isobaric_potential = Quantity(-74150 * units.joule / units.mole)
    temperature = Quantity(2000 * units.kelvin)

    return Args(standart_change_isobaric_potential=standart_change_isobaric_potential,
        temperature=temperature)


def test_basic_equilibrium_constant(test_args: Args) -> None:
    result = constant_law.calculate_equilibrium_constant(test_args.standart_change_isobaric_potential, test_args.temperature)
    assert_equal(result, 86.409)


def test_bad_standart_change_isobaric_potential(test_args: Args) -> None:
    standart_change_isobaric_potential = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        constant_law.calculate_equilibrium_constant(standart_change_isobaric_potential, test_args.temperature)
    with raises(TypeError):
        constant_law.calculate_equilibrium_constant(100, test_args.temperature)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        constant_law.calculate_equilibrium_constant(test_args.standart_change_isobaric_potential, temperature)
    with raises(TypeError):
        constant_law.calculate_equilibrium_constant(test_args.standart_change_isobaric_potential, 100)
