from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as movement_law

# Description
## The man starts running with initial speed 10m/s and is getting tired while running (decreases his velocity) with -0.1 m/s/s. How far this man can run away in 10 seconds?

Args = namedtuple("Args", ["V0", "a", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    V0 = Quantity(10 * units.meter / units.second)
    a = Quantity(-0.1 * units.meter / units.second**2)
    t = Quantity(10 * units.second)
    return Args(V0=V0, a=a, t=t)


def test_basic_distance(test_args: Args) -> None:
    result = movement_law.calculate_distance(test_args.V0, test_args.a, test_args.t)
    assert_equal(result, 95 * units.meter)


def test_bad_velocity(test_args: Args) -> None:
    Vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(Vb, test_args.a, test_args.t)
    with raises(TypeError):
        movement_law.calculate_distance(100, test_args.a, test_args.t)


def test_bad_acceleration(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(test_args.V0, ab, test_args.t)
    with raises(TypeError):
        movement_law.calculate_distance(test_args.V0, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(test_args.V0, test_args.a, tb)
    with raises(TypeError):
        movement_law.calculate_distance(test_args.V0, test_args.a, 100)
