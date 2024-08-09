from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    SI,
)
from symplyphysics.laws.thermodynamics import total_energy_transfer_is_zero_in_isolated_system as heat_balance_law

Args = namedtuple("Args", ["q1", "q2", "q3"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q1 = Quantity(3 * units.joules)
    q2 = Quantity(-5 * units.joules)
    q3 = Quantity(2 * units.joules)
    return Args(q1=q1, q2=q2, q3=q3)


def test_basic_amount_heat_balance(test_args: Args) -> None:
    result = heat_balance_law.calculate_amount_energy([test_args.q1, test_args.q2])
    assert_equal(result, 2 * units.joule)


def test_three_array_elements(test_args: Args) -> None:
    result = heat_balance_law.calculate_amount_energy([test_args.q1, test_args.q2, test_args.q3])
    assert_equal(result, 0 * units.joule)


def test_array_empty() -> None:
    result = heat_balance_law.calculate_amount_energy([])
    assert SI.get_dimension_system().is_dimensionless(result.dimension)
    assert_equal(result, 0)


def test_array_with_bad_element(test_args: Args) -> None:
    qb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        heat_balance_law.calculate_amount_energy([test_args.q1, qb])
    with raises(TypeError):
        heat_balance_law.calculate_amount_energy([test_args.q1, 100])
    with raises(errors.UnitsError):
        heat_balance_law.calculate_amount_energy([qb, test_args.q2])
    with raises(TypeError):
        heat_balance_law.calculate_amount_energy([100, test_args.q2])
    with raises(errors.UnitsError):
        heat_balance_law.calculate_amount_energy([qb, qb])
    with raises(TypeError):
        heat_balance_law.calculate_amount_energy([100, 100])
    with raises(TypeError):
        heat_balance_law.calculate_amount_energy(test_args.q1)
