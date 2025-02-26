from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.relativistic import lorentz_transformation_of_time as time_law

Args = namedtuple("Args", "t x v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(1 * units.second)
    x = Quantity(100 * units.kilometer)
    v = Quantity(0.5 * units.speed_of_light)
    return Args(t=t, x=x, v=v)


def test_law(test_args: Args) -> None:
    result = time_law.calculate_time_in_proper_frame(test_args.t, test_args.x, test_args.v)
    assert_equal(result, 1.15 * units.second, relative_tolerance=4e-3)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        time_law.calculate_time_in_proper_frame(tb, test_args.x, test_args.v)
    with raises(TypeError):
        time_law.calculate_time_in_proper_frame(100, test_args.x, test_args.v)


def test_bad_position(test_args: Args) -> None:
    xb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        time_law.calculate_time_in_proper_frame(test_args.t, xb, test_args.v)
    with raises(TypeError):
        time_law.calculate_time_in_proper_frame(test_args.t, 100, test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        time_law.calculate_time_in_proper_frame(test_args.t, test_args.x, vb)
    with raises(TypeError):
        time_law.calculate_time_in_proper_frame(test_args.t, test_args.x, 100)
