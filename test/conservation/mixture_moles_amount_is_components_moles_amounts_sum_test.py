from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.conservation import mixture_moles_amount_is_components_moles_amounts_sum as sum_of_moles_law

Args = namedtuple("Args", ["nu1", "nu2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    nu1 = Quantity(10 * units.moles)
    nu2 = Quantity(20 * units.moles)
    return Args(nu1=nu1, nu2=nu2)


def test_basic_amount_of_moles(test_args: Args) -> None:
    result = sum_of_moles_law.calculate_moles_count_of_mixture([test_args.nu1, test_args.nu2])
    assert_equal(result, 30 * units.moles)


def test_three_amounts_of_moles_array(test_args: Args) -> None:
    m3 = Quantity(5 * units.moles)
    result = sum_of_moles_law.calculate_moles_count_of_mixture([test_args.nu1, test_args.nu2, m3])
    assert_equal(result, 35 * units.moles)


def test_bad_amount_of_moles(test_args: Args) -> None:
    nu_b = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        sum_of_moles_law.calculate_moles_count_of_mixture([nu_b, test_args.nu2])
    with raises(TypeError):
        sum_of_moles_law.calculate_moles_count_of_mixture([100, test_args.nu2])
    with raises(errors.UnitsError):
        sum_of_moles_law.calculate_moles_count_of_mixture([test_args.nu1, nu_b])
    with raises(TypeError):
        sum_of_moles_law.calculate_moles_count_of_mixture([test_args.nu1, 100])
    with raises(errors.UnitsError):
        sum_of_moles_law.calculate_moles_count_of_mixture([nu_b, nu_b])
    with raises(TypeError):
        sum_of_moles_law.calculate_moles_count_of_mixture([100, 100])
    with raises(TypeError):
        sum_of_moles_law.calculate_moles_count_of_mixture(test_args.nu1)
