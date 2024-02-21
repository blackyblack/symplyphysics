from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import damped_harmonic_oscillator_equation as damped_eqn

# Description
## A damped oscillator is critically damped, with the damping constant z = 2.5, and an
## undamped angular frequency w = 5 Hz. Initial position of the oscillator is 0, and
## its initial speed is 1 m/s. Then at t = 0.5 s, the oscillator's position is 0.041 m.

Args = namedtuple("Args", "x0 v0 w z t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    x0 = Quantity(0.0 * units.meter)
    v0 = Quantity(1.0 * units.meter / units.second)
    w = Quantity(5.0 * units.hertz)
    z = 2.5
    t = Quantity(0.5 * units.second)
    return Args(x0=x0, v0=v0, w=w, z=z, t=t)


def test_critically_damped(test_args: Args) -> None:
    result = damped_eqn.calculate_displacement(
        test_args.x0, test_args.v0, test_args.w, test_args.z, test_args.t
    )
    assert_equal(result, 0.041 * units.meter, tolerance=2e-3)


def test_bad_initial_position(test_args: Args) -> None:
    xb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        damped_eqn.calculate_displacement(
            xb, test_args.v0, test_args.w, test_args.z, test_args.t
        )
    with raises(TypeError):
        damped_eqn.calculate_displacement(
            100, test_args.v0, test_args.w, test_args.z, test_args.t
        )


def test_bad_initial_velocity(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        damped_eqn.calculate_displacement(
            test_args.x0, vb, test_args.w, test_args.z, test_args.t
        )
    with raises(TypeError):
        damped_eqn.calculate_displacement(
            test_args.x0, 100, test_args.w, test_args.z, test_args.t
        )


def test_bad_undamped_angular_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        damped_eqn.calculate_displacement(
            test_args.x0, test_args.v0, wb, test_args.z, test_args.t
        )
    with raises(TypeError):
        damped_eqn.calculate_displacement(
            test_args.x0, test_args.v0, 100, test_args.z, test_args.t
        )

def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        damped_eqn.calculate_displacement(
            test_args.x0, test_args.v0, test_args.w, test_args.z, tb
        )
    with raises(TypeError):
        damped_eqn.calculate_displacement(
            test_args.x0, test_args.v0, test_args.w, test_args.z, 100
        )
