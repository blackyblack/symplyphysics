from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.electricity import power_from_energy_time as power_def


# How much power did the heater use if it is known that it gave off 20 kilojoules
# of energy in 35 seconds? Consider that all energy consumed equals energy given up.
@fixture(name="test_args")
def test_args_fixture():
    Q = Quantity(20 * prefixes.kilo * units.joule)
    t = Quantity(35 * units.second)
    Args = namedtuple("Args", ["Q", "t"])
    return Args(Q=Q, t=t)


def test_basic_power(test_args):
    result = power_def.calculate_power(test_args.Q, test_args.t)
    assert_equal(result, 571 * units.watt)


def test_bad_energy(test_args):
    Qb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        power_def.calculate_power(Qb, test_args.t)
    with raises(TypeError):
        power_def.calculate_power(100, test_args.t)


def test_bad_time(test_args):
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        power_def.calculate_power(test_args.Q, tb)
    with raises(TypeError):
        power_def.calculate_power(test_args.Q, 100)
