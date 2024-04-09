from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity import second_escape_velocity

# Second cosmic velocity for Earth near surface is 11,2 km/s

Args = namedtuple("Args", ["m", "h", "r"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(5.97e24 * units.kilograms)
    r = Quantity(6_371 * units.kilometers)
    h = Quantity(0 * units.meters)
    return Args(m=m, h=h, r=r)


def test_basic_velocity(test_args: Args) -> None:
    result = second_escape_velocity.calculate_velocity(test_args.m, test_args.r, test_args.h)
    assert_equal(result, 11.18 * units.kilometer / units.second)


def test_bad_mass(test_args: Args) -> None:
    bm = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        second_escape_velocity.calculate_velocity(bm, test_args.r, test_args.h)
    with raises(TypeError):
        second_escape_velocity.calculate_velocity(100, test_args.r, test_args.h)


def test_bad_radius(test_args: Args) -> None:
    br = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        second_escape_velocity.calculate_velocity(test_args.m, br, test_args.h)
    with raises(TypeError):
        second_escape_velocity.calculate_velocity(test_args.m, 100, test_args.h)


def test_bad_height(test_args: Args) -> None:
    bh = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        second_escape_velocity.calculate_velocity(test_args.m, test_args.r, bh)
    with raises(TypeError):
        second_escape_velocity.calculate_velocity(test_args.m, test_args.r, 100)
