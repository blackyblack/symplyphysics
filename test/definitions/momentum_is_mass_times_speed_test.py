from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import momentum_is_mass_times_speed as momentum_def

Args = namedtuple("Args", ["m", "v"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(1 * units.kilogram)
    v = Quantity(5 * units.meter / units.second)
    return Args(m=m, v=v)


def test_basic_momentum(test_args: Args) -> None:
    result = momentum_def.calculate_momentum(test_args.m, test_args.v)
    assert_equal(result, 5 * units.kilogram * units.meter / units.second)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(mb, test_args.v)
    with raises(TypeError):
        momentum_def.calculate_momentum(100, test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(test_args.m, vb)
    with raises(TypeError):
        momentum_def.calculate_momentum(test_args.m, 100)
