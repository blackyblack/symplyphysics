from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_approx, Quantity, units, errors
from symplyphysics.laws.hydro import mach_number

# Example from: https://www.omnicalculator.com/physics/mach-number


@fixture(name="test_args")
def test_args_fixture():
    velocity = Quantity(120 * units.kilometer / units.hour)
    speed_of_sound = Quantity(343.2 * units.meter / units.second)
    Args = namedtuple("Args", ["velocity", "speed_of_sound"])
    return Args(velocity=velocity, speed_of_sound=speed_of_sound)


def test_basic_mach_number(test_args):
    result = mach_number.calculate_mach_number(test_args.velocity, test_args.speed_of_sound)
    assert_approx(result, 0.09712)


def test_bad_velocity(test_args):
    bv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mach_number.calculate_mach_number(bv, test_args.speed_of_sound)
    with raises(TypeError):
        mach_number.calculate_mach_number(0, test_args.speed_of_sound)


def test_bad_speed_of_sound(test_args):
    bss = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mach_number.calculate_mach_number(test_args.velocity, bss)
    with raises(TypeError):
        mach_number.calculate_mach_number(test_args.velocity, 0)
