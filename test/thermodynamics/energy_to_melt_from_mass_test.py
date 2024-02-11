from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.thermodynamics import energy_to_melt_from_mass as amount_energy

# https://easyfizika.ru/zadachi/termodinamika/vannu-emkostyu-100-litrov-neobhodimo-zapolnit-vodoj-imeyushhej-temperaturu-30-c/
# If mass of ice equal 29.7 kg and specific heat of melting equal 330 KJ/kg,
# the melting energy should be equal to 9.801 MJ

Args = namedtuple("Args", ["k_lambda", "m"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    k_lambda = Quantity(330 * prefixes.kilo * units.joule / units.kilogram)
    m = Quantity(29.7 * units.kilogram)
    return Args(k_lambda=k_lambda, m=m)


def test_basic_amount(test_args: Args) -> None:
    result = amount_energy.calculate_amount_energy(test_args.k_lambda, test_args.m)
    assert_equal(result, 9.801 * prefixes.mega * units.joules)


def test_bad_specific_heat_melting(test_args: Args) -> None:
    k_lambdab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(k_lambdab, test_args.m)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(100, test_args.m)


def test_bad_body_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.k_lambda, mb)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.k_lambda, 100)
