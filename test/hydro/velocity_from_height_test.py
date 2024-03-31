from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.hydro import velocity_from_height as torricellis_formula

# Description
## In the bottom of water tank there is a small hole. Height of liquid in tank is 3m. Velocity of jet according to Torricelli's formula should be
## 7.67m/s

Args = namedtuple("Args", ["h"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    h = Quantity(3 * units.meter)
    return Args(h=h)


def test_basic_velocity(test_args: Args) -> None:
    result = torricellis_formula.calculate_velocity(test_args.h)
    assert_equal(result, 7.67 * units.meter / units.second)


def test_bad_height() -> None:
    hb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        torricellis_formula.calculate_velocity(hb)
    with raises(TypeError):
        torricellis_formula.calculate_velocity(100)
