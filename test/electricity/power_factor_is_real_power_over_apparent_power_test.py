from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.laws.electricity import power_factor_is_real_power_over_apparent_power as power_factor_law

# Description
## Assert we have a device wich consumes 10 Watt of power and makes 3 Watt of work. Power factor of this consumer should be 3/10.

Args = namedtuple("Args", ["P", "S"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    P = Quantity(3 * units.watt)
    S = Quantity(10 * units.watt)
    return Args(P=P, S=S)


def test_basic_power_factor(test_args: Args) -> None:
    result = power_factor_law.calculate_power_factor(test_args.P, test_args.S)
    assert_equal(result, 0.3)


def test_bad_power(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_factor_law.calculate_power_factor(pb, test_args.S)
        power_factor_law.calculate_power_factor(test_args.P, pb)
    with raises(TypeError):
        power_factor_law.calculate_power_factor(100, test_args.S)
        power_factor_law.calculate_power_factor(test_args.P, 100)
