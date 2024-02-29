from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic.damped_oscillations import critical_damping

# Description
## An damped oscillator is in the critically damped state. Its angular frequency
## when it is undamped is 1 Hz. Initially it was at position x0 = 0.2 m and had
## no initial velocity. Then it was at position x = 8.12 cm at time t = 2 s.

Args = namedtuple("Args", "x0 v0 omega t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    x0 = Quantity(0.2 * units.meter)
    v0 = Quantity(0 * units.meter / units.second)
    omega = Quantity(1.0 * units.hertz)
    t = Quantity(2.0 * units.second)
    return Args(x0=x0, v0=v0, omega=omega, t=t)


def test_law(test_args: Args) -> None:
    result = critical_damping.calculate_displacement(test_args.x0, test_args.v0, test_args.omega,
        test_args.t)
    assert_equal(result, 8.12 * units.centimeter)


def test_bad_initial_position(test_args: Args) -> None:
    xb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        critical_damping.calculate_displacement(xb, test_args.v0, test_args.omega, test_args.t)
    with raises(TypeError):
        critical_damping.calculate_displacement(100, test_args.v0, test_args.omega, test_args.t)


def test_bad_initial_velocity(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        critical_damping.calculate_displacement(test_args.x0, vb, test_args.omega, test_args.t)
    with raises(TypeError):
        critical_damping.calculate_displacement(test_args.x0, 100, test_args.omega, test_args.t)


def test_bad_angular_frequency(test_args: Args) -> None:
    omega_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        critical_damping.calculate_displacement(test_args.x0, test_args.v0, omega_bad, test_args.t)
    with raises(TypeError):
        critical_damping.calculate_displacement(test_args.x0, test_args.v0, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        critical_damping.calculate_displacement(test_args.x0, test_args.v0, test_args.omega, tb)
    with raises(TypeError):
        critical_damping.calculate_displacement(test_args.x0, test_args.v0, test_args.omega, 100)
