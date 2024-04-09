from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import relationship_of_periods_of_planets_with_mass_and_semi_axes_of_orbits as period_law

# Description
## The semimajor axis of the Earth is 1 astronomical unit, and the semimajor axis of the Mars is 1.524 astronomical unit.
## The rotation period of the Mars is 1.88 year. Then the rotation period of the Earth will be 1 year.

Args = namedtuple("Args",
    ["second_period", "first_semi_axis", "second_semi_axis", "first_mass", "second_mass"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    second_period = Quantity(1.88 * units.year)
    first_semi_axis = Quantity(1 * units.astronomical_unit)
    second_semi_axis = Quantity(1.524 * units.astronomical_unit)
    first_mass = Quantity(5.9742e24 * units.kilogram)
    second_mass = Quantity(6.39e23 * units.kilogram)

    return Args(second_period=second_period,
        first_semi_axis=first_semi_axis,
        second_semi_axis=second_semi_axis,
        first_mass=first_mass,
        second_mass=second_mass)


def test_basic_first_period(test_args: Args) -> None:
    result = period_law.calculate_first_period(test_args.second_period, test_args.first_semi_axis,
        test_args.second_semi_axis, test_args.first_mass, test_args.second_mass)
    assert_equal(result, 1 * units.year)


def test_bad_second_period(test_args: Args) -> None:
    second_period = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        period_law.calculate_first_period(second_period, test_args.first_semi_axis,
            test_args.second_semi_axis, test_args.first_mass, test_args.second_mass)
    with raises(TypeError):
        period_law.calculate_first_period(100, test_args.first_semi_axis,
            test_args.second_semi_axis, test_args.first_mass, test_args.second_mass)


def test_bad_semi_axis(test_args: Args) -> None:
    bad_semi_axis = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        period_law.calculate_first_period(test_args.second_period, bad_semi_axis,
            test_args.second_semi_axis, test_args.first_mass, test_args.second_mass)
    with raises(TypeError):
        period_law.calculate_first_period(test_args.second_period, 100, test_args.second_semi_axis,
            test_args.first_mass, test_args.second_mass)
    with raises(errors.UnitsError):
        period_law.calculate_first_period(test_args.second_period, test_args.first_semi_axis,
            bad_semi_axis, test_args.first_mass, test_args.second_mass)
    with raises(TypeError):
        period_law.calculate_first_period(test_args.second_period, test_args.first_semi_axis, 100,
            test_args.first_mass, test_args.second_mass)


def test_bad_mass(test_args: Args) -> None:
    bad_mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        period_law.calculate_first_period(test_args.second_period, test_args.first_semi_axis,
            test_args.second_semi_axis, bad_mass, test_args.second_mass)
    with raises(TypeError):
        period_law.calculate_first_period(test_args.second_period, test_args.first_semi_axis,
            test_args.second_semi_axis, 100, test_args.second_mass)
    with raises(errors.UnitsError):
        period_law.calculate_first_period(test_args.second_period, test_args.first_semi_axis,
            test_args.second_semi_axis, test_args.first_mass, bad_mass)
    with raises(TypeError):
        period_law.calculate_first_period(test_args.second_period, test_args.first_semi_axis,
            test_args.second_semi_axis, test_args.first_mass, 100)
