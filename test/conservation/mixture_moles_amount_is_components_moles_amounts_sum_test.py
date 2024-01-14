from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.core import errors
from symplyphysics.laws.conservation import mixture_moles_amount_is_components_moles_amounts_sum as sum_of_moles_law


@fixture(name="test_args")
def test_args_fixture():
    nu1 = Quantity(10 * units.moles)
    nu2 = Quantity(20 * units.moles)
    Args = namedtuple("Args", ["nu1", "nu2"])
    return Args(nu1=nu1, nu2=nu2)


def test_basic_amount_of_moles(test_args):
    result = sum_of_moles_law.calculate_moles_count_of_mixture([test_args.nu1, test_args.nu2])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.amount_of_substance)
    result_amount_of_moles = convert_to(result, units.moles).evalf(3)
    assert result_amount_of_moles == approx(30, 0.001)


def test_three_masses_array(test_args):
    m3 = Quantity(5 * units.moles)
    result = sum_of_moles_law.calculate_moles_count_of_mixture([test_args.nu1, test_args.nu2, m3])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.amount_of_substance)
    result_amount_of_moles = convert_to(result, units.moles).evalf(3)
    assert result_amount_of_moles == approx(35, 0.001)


def test_bad_amount_of_moles(test_args):
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
