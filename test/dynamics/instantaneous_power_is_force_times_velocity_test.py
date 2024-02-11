from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import (
    instantaneous_power_is_force_times_velocity as power_law,)

# Description
## A particle with velocity v = 2.45 m/s moves due to a force F = 2.24 N, the angle between
## the two vectors being 0.752 rad. The instantaneous power of the force amounts to 4.01 W.

Args = namedtuple("Args", "v f phi")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v = Quantity(2.45 * units.meter / units.second)
    f = Quantity(2.24 * units.newton)
    phi = Quantity(0.752 * units.radian)
    return Args(v=v, f=f, phi=phi)


def test_law(test_args: Args) -> None:
    result = power_law.calculate_power(test_args.f, test_args.v, test_args.phi)
    assert_equal(result, 4.01 * units.watt)


def test_bad_force(test_args: Args) -> None:
    fb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        power_law.calculate_power(fb, test_args.v, test_args.phi)
    with raises(TypeError):
        power_law.calculate_power(100, test_args.v, test_args.phi)


def test_bad_speed(test_args: Args) -> None:
    vb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        power_law.calculate_power(test_args.f, vb, test_args.phi)
    with raises(TypeError):
        power_law.calculate_power(test_args.f, 100, test_args.phi)


def test_bad_angle(test_args: Args) -> None:
    phib = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        power_law.calculate_power(test_args.f, test_args.v, phib)
