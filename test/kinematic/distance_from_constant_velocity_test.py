from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematics import distance_from_constant_velocity as movement_law

# Description
## If an object is traveling at a constant velocity of 12m/s, then calculate the distance covered by the object after 1 minute.

Args = namedtuple("Args", ["x0", "v", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    x0 = Quantity(0 * units.meter)
    v = Quantity(12 * units.meter / units.second)
    t = Quantity(60 * units.second)
    return Args(x0=x0, v=v, t=t)


def test_basic_distance(test_args: Args) -> None:
    result = movement_law.calculate_distance(test_args.x0, test_args.v, test_args.t)
    assert_equal(result, 720 * units.meter)


def test_bad_distance(test_args: Args) -> None:
    xb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(xb, test_args.v, test_args.t)
    with raises(TypeError):
        movement_law.calculate_distance(100, test_args.v, test_args.t)


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(test_args.x0, vb, test_args.t)
    with raises(TypeError):
        movement_law.calculate_distance(test_args.x0, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(test_args.x0, test_args.v, tb)
    with raises(TypeError):
        movement_law.calculate_distance(test_args.x0, test_args.v, 100)
