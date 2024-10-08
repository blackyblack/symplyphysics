from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import maximum_height_from_initial_speed as maximum_height_law

# Description
## If a body is thrown vertically upwards at a speed of 10 m/s, it will rise to a height of  (10*10) /2*9,8 = 5.1 meters

Args = namedtuple("Args", ["v"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v = Quantity(10 * units.meter / units.second)
    return Args(v=v)


def test_basic_height(test_args: Args) -> None:
    result = maximum_height_law.calculate_maximum_height(test_args.v)
    assert_equal(result, 5.1 * units.meter)


def test_bad_velocity() -> None:
    vb = Quantity(1 * units.coulomb)

    with raises(errors.UnitsError):
        maximum_height_law.calculate_maximum_height(vb)

    with raises(TypeError):
        maximum_height_law.calculate_maximum_height(100)
