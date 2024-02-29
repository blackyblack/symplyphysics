from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic import displacement_in_simple_harmonic_motion as harmonic_law

# Description
## An force of oscillating magnitude is acting on a particle. The amplitude of its oscillatons is 1.0 N,
## their angular frequency is 500 rad/s and the phase lag is zero. At t = 0.1 s, the magnitude
## of the force amounts to 0.965 N.

Args = namedtuple("Args", "f_m w t phi")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f_m = Quantity(1.0 * units.newton)
    w = Quantity(500.0 * units.radian / units.second)
    phi = 0.0
    t = Quantity(0.1 * units.second)
    return Args(f_m=f_m, w=w, t=t, phi=phi)


def test_law(test_args: Args) -> None:
    result = harmonic_law.calculate_displacement(test_args.f_m, test_args.w, test_args.phi, test_args.t)
    assert_equal(result, 0.965 * units.newton)


def test_bad_angular_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        harmonic_law.calculate_displacement(test_args.f_m, wb, test_args.phi, test_args.t)
    with raises(TypeError):
        harmonic_law.calculate_displacement(test_args.f_m, 100, test_args.phi, test_args.t)


def test_bad_phase_lag(test_args: Args) -> None:
    phi_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        harmonic_law.calculate_displacement(test_args.f_m, test_args.w, phi_bad, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        harmonic_law.calculate_displacement(test_args.f_m, test_args.w, test_args.phi, tb)
    with raises(TypeError):
        harmonic_law.calculate_displacement(test_args.f_m, test_args.w, test_args.phi, 100)
