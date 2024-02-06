from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity import escape_velocity

# First cosmic velocity for Earth near surface is 7,91 km/s


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(5.97e24 * units.kilograms)
    r = Quantity(6_371 * units.kilometers)
    h = Quantity(0 * units.meters)
    Args = namedtuple("Args", ["m", "h", "r"])
    return Args(m=m, h=h, r=r)


def test_basic_velocity(test_args):
    result = escape_velocity.calculate_velocity(test_args.m, test_args.r, test_args.h)
    assert_equal(result, 7910 * units.meter / units.second)


def test_bad_mass(test_args):
    bm = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        escape_velocity.calculate_velocity(bm, test_args.r, test_args.h)
    with raises(TypeError):
        escape_velocity.calculate_velocity(100, test_args.r, test_args.h)


def test_bad_radius(test_args):
    br = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        escape_velocity.calculate_velocity(test_args.m, br, test_args.h)
    with raises(TypeError):
        escape_velocity.calculate_velocity(test_args.m, 100, test_args.h)


def test_bad_height(test_args):
    bh = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        escape_velocity.calculate_velocity(test_args.m, test_args.r, bh)
    with raises(TypeError):
        escape_velocity.calculate_velocity(test_args.m, test_args.r, 100)
