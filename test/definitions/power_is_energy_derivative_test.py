from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import power_is_energy_derivative as power_def

Args = namedtuple("Args", ["Q0", "Q1", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    Q0 = Quantity(0 * units.joule)
    Q1 = Quantity(20 * units.joule)
    t = Quantity(5 * units.second)
    return Args(Q0=Q0, Q1=Q1, t=t)


def test_basic_power(test_args: Args) -> None:
    result = power_def.calculate_power(test_args.Q0, test_args.Q1, test_args.t)
    assert_equal(result, 4 * units.watt)


def test_power_with_bad_energy(test_args: Args) -> None:
    Qb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        power_def.calculate_power(Qb, test_args.Q1, test_args.t)
    with raises(errors.UnitsError):
        power_def.calculate_power(test_args.Q0, Qb, test_args.t)
    with raises(TypeError):
        power_def.calculate_power(100, test_args.Q1, test_args.t)
    with raises(TypeError):
        power_def.calculate_power(test_args.Q0, 100, test_args.t)


def test_power_with_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        power_def.calculate_power(test_args.Q0, test_args.Q1, tb)
    with raises(TypeError):
        power_def.calculate_power(test_args.Q0, test_args.Q1, 100)
