from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import transient_attenuation_of_separate_loop_coupler_with_cascade_connection as attenuation_law

# Description
## The required transient attenuation of the cascade is -9 decibels, number of couplers is 6.
## Then the transient attenuation of one coupler is equal to -1.7 decibels.

Args = namedtuple("Args", ["attenuation_of_cascade", "number_of_couplers"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    attenuation_of_cascade = -9
    number_of_couplers = 6

    return Args(attenuation_of_cascade=attenuation_of_cascade,
        number_of_couplers=number_of_couplers)


def test_basic_attenuation(test_args: Args) -> None:
    result = attenuation_law.calculate_attenuation_of_coupler(test_args.attenuation_of_cascade,
        test_args.number_of_couplers)
    assert_equal(result, -1.7)


def test_bad_attenuation_of_cascade(test_args: Args) -> None:
    attenuation_of_cascade = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        attenuation_law.calculate_attenuation_of_coupler(attenuation_of_cascade, test_args.number_of_couplers)


def test_bad_number_of_couplers(test_args: Args) -> None:
    number_of_couplers = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        attenuation_law.calculate_attenuation_of_coupler(test_args.attenuation_of_cascade, number_of_couplers)
