from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import relationship_of_periods_of_planets_with_semi_axes_of_orbits as period_law

# Description
## Let the semimajor axis of the first planet be 784122156 kilometer, and the semimajor axis of the second planet be 149599300 kilometer.
## The rotation period of the second planet is 1 year. Then the rotation period of the first planet will be 12 years.
## https://infourok.ru/astronomiya-reshenie-zadach-po-teme-zakony-keplera-4414880.html

Args = namedtuple("Args", ["second_period", "first_semi_axis", "second_semi_axis"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    second_period = Quantity(1 * units.year)
    first_semi_axis = Quantity(784122156 * units.kilometer)
    second_semi_axis = Quantity(149599300 * units.kilometer)

    return Args(second_period=second_period,
        first_semi_axis=first_semi_axis,
        second_semi_axis=second_semi_axis)


def test_basic_first_period(test_args: Args) -> None:
    result = period_law.calculate_first_period(test_args.second_period,
        test_args.first_semi_axis, test_args.second_semi_axis)
    assert_equal(result, 12 * units.year)


def test_bad_second_period(test_args: Args) -> None:
    second_period = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        period_law.calculate_first_period(second_period, test_args.first_semi_axis,
            test_args.second_semi_axis)
    with raises(TypeError):
        period_law.calculate_first_period(100, test_args.first_semi_axis, test_args.second_semi_axis)


def test_bad_bad_semi_axis(test_args: Args) -> None:
    bad_semi_axis = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        period_law.calculate_first_period(test_args.second_period, bad_semi_axis,
            test_args.second_semi_axis)
    with raises(TypeError):
        period_law.calculate_first_period(test_args.second_period, 100,
            test_args.second_semi_axis)
    with raises(errors.UnitsError):
        period_law.calculate_first_period(test_args.second_period,
            test_args.first_semi_axis, bad_semi_axis)
    with raises(TypeError):
        period_law.calculate_first_period(test_args.second_period,
            test_args.first_semi_axis, 100)
