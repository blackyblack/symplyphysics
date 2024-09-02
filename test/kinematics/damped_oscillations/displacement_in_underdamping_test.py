from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematics.damped_oscillations import displacement_in_underdamping

# Description
## An oscillator is underdamped, its undamped angular frequency w = 5 Hz, exponential
## decay constant lambda = 0.3 1/s, and the phase lag phi = 0.1 rad. The initial
## position of the oscillator is 0.5 m. At time t = 2 s, its position is x = -0.215 m.

Args = namedtuple("Args", "x_init lambda_ w phi t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    x_init = Quantity(0.5 * units.meter)
    lambda_ = Quantity(0.3 / units.second)
    w = Quantity(5.0 * units.hertz)
    phi = Quantity(0.1 * units.radian)
    t = Quantity(2.0 * units.second)
    return Args(x_init=x_init, lambda_=lambda_, w=w, phi=phi, t=t)


def test_law(test_args: Args) -> None:
    result = displacement_in_underdamping.calculate_displacement(test_args.x_init,
        test_args.lambda_, test_args.w, test_args.phi, test_args.t)
    assert_equal(result, -0.215 * units.meter, tolerance=2e-3)


def test_bad_initial_position(test_args: Args) -> None:
    xb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        displacement_in_underdamping.calculate_displacement(xb, test_args.lambda_, test_args.w,
            test_args.phi, test_args.t)
    with raises(TypeError):
        displacement_in_underdamping.calculate_displacement(100, test_args.lambda_, test_args.w,
            test_args.phi, test_args.t)


def test_bad_decay_constant(test_args: Args) -> None:
    lambda_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        displacement_in_underdamping.calculate_displacement(test_args.x_init, lambda_bad,
            test_args.w, test_args.phi, test_args.t)
    with raises(TypeError):
        displacement_in_underdamping.calculate_displacement(test_args.x_init, 100, test_args.w,
            test_args.phi, test_args.t)


def test_bad_undamped_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        displacement_in_underdamping.calculate_displacement(test_args.x_init, test_args.lambda_, wb,
            test_args.phi, test_args.t)
    with raises(TypeError):
        displacement_in_underdamping.calculate_displacement(test_args.x_init, test_args.lambda_,
            100, test_args.phi, test_args.t)


def test_bad_phase_lag(test_args: Args) -> None:
    phi_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        displacement_in_underdamping.calculate_displacement(test_args.x_init, test_args.lambda_,
            test_args.w, phi_bad, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        displacement_in_underdamping.calculate_displacement(test_args.x_init, test_args.lambda_,
            test_args.w, test_args.phi, tb)
    with raises(TypeError):
        displacement_in_underdamping.calculate_displacement(test_args.x_init, test_args.lambda_,
            test_args.w, test_args.phi, 100)
