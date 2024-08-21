from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, Quantity, units, errors
from symplyphysics.laws.hydro import mach_number_is_flow_speed_over_speed_of_sound

# Example from: https://www.omnicalculator.com/physics/mach-number

Args = namedtuple("Args", ["velocity", "speed_of_sound"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    velocity = Quantity(120 * units.kilometer / units.hour)
    speed_of_sound = Quantity(343.2 * units.meter / units.second)
    return Args(velocity=velocity, speed_of_sound=speed_of_sound)


def test_basic_mach_number(test_args: Args) -> None:
    result = mach_number_is_flow_speed_over_speed_of_sound.calculate_mach_number(test_args.velocity, test_args.speed_of_sound)
    assert_equal(result, 0.09712)


def test_bad_velocity(test_args: Args) -> None:
    bv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mach_number_is_flow_speed_over_speed_of_sound.calculate_mach_number(bv, test_args.speed_of_sound)
    with raises(TypeError):
        mach_number_is_flow_speed_over_speed_of_sound.calculate_mach_number(0, test_args.speed_of_sound)


def test_bad_speed_of_sound(test_args: Args) -> None:
    bss = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mach_number_is_flow_speed_over_speed_of_sound.calculate_mach_number(test_args.velocity, bss)
    with raises(TypeError):
        mach_number_is_flow_speed_over_speed_of_sound.calculate_mach_number(test_args.velocity, 0)
