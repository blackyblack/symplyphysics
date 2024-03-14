from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import mach_number as mach_number_def

# Description
## The speed of sound in the air at normal conditions is 343 m/s. An object is moving
## in the air at a speed of 1000 m/s. The Mach number of the moving object is 2.9

Args = namedtuple("Args", "u c")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    u = Quantity(1000 * units.meter / units.second)
    c = Quantity(343 * units.meter / units.second)
    return Args(u=u, c=c)


def test_law(test_args: Args) -> None:
    result = mach_number_def.calculate_mach_number(test_args.u, test_args.c)
    assert_equal(result, 2.9, tolerance=6e-3)


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        mach_number_def.calculate_mach_number(vb, test_args.c)
    with raises(TypeError):
        mach_number_def.calculate_mach_number(100, test_args.c)
    with raises(errors.UnitsError):
        mach_number_def.calculate_mach_number(test_args.u, vb)
    with raises(TypeError):
        mach_number_def.calculate_mach_number(test_args.u, 100)
