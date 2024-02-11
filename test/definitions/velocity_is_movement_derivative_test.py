from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import velocity_is_movement_derivative as velocity_def

# Description
## Assume object starts moving with some velocity. After 5 seconds it's distance is 80 meters.
## Velocity should be 16 m/s.

Args = namedtuple("Args", ["S0", "S1", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    S0 = Quantity(0 * units.meter)
    S1 = Quantity(80 * units.meters)
    t = Quantity(5 * units.second)
    return Args(S0=S0, S1=S1, t=t)


def test_basic_velocity(test_args: Args) -> None:
    result = velocity_def.calculate_velocity(test_args.S0, test_args.S1, test_args.t)
    assert_equal(result, 16 * units.meter / units.second)


def test_velocity_with_bad_distance(test_args: Args) -> None:
    Sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_def.calculate_velocity(Sb, test_args.S1, test_args.t)
    with raises(errors.UnitsError):
        velocity_def.calculate_velocity(test_args.S0, Sb, test_args.t)
    with raises(TypeError):
        velocity_def.calculate_velocity(100, test_args.S1, test_args.t)
    with raises(TypeError):
        velocity_def.calculate_velocity(test_args.S0, 100, test_args.t)


def test_velocity_with_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_def.calculate_velocity(test_args.S0, test_args.S1, tb)
    with raises(TypeError):
        velocity_def.calculate_velocity(test_args.S0, test_args.S1, 100)
