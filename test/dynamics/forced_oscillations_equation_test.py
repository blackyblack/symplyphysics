from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import forced_oscillations_equation as forced_eqn

# Description
## An oscillating external force is acting on an oscillator of mass m = 1 kg and natural
## angular frequency w0 = 1 rad/s. The amplitude of the driving force is f = 5 N, its angular
## frequency is 0.95 rad/s and its phase lag is pi/3. At t = 15 s, the displacement of the
## oscillator amounts to -6.46 cm. Initially the oscillator was 5 cm displaced and had a speed
## of -1 cm/s.

Args = namedtuple("Args", "q0 v0 m w0 f w phi t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q0 = Quantity(5.0 * units.centimeter)
    v0 = Quantity(-1.0 * units.centimeter / units.second)
    m = Quantity(1.0 * units.kilogram)
    w0 = Quantity(1.0 * units.radian / units.second)
    f = Quantity(1.0 * units.newton)
    w = Quantity(0.95 * units.radian / units.second)
    phi = Quantity(pi / 3 * units.radian)
    t = Quantity(15.0 * units.second)
    return Args(q0=q0, v0=v0, m=m, w0=w0, f=f, w=w, phi=phi, t=t)


def test_law(test_args: Args) -> None:
    result = forced_eqn.calculate_displacement(
        test_args.q0, test_args.v0, test_args.m, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
    )
    assert_equal(result, -6.46 * units.centimeter)


def test_bad_position(test_args: Args) -> None:
    qb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_eqn.calculate_displacement(
            qb, test_args.v0, test_args.m, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
        )
    with raises(TypeError):
        forced_eqn.calculate_displacement(
            100, test_args.v0, test_args.m, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
        )


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_eqn.calculate_displacement(
            test_args.q0, vb, test_args.m, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
        )
    with raises(TypeError):
        forced_eqn.calculate_displacement(
            test_args.q0, 100, test_args.m, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
        )


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_eqn.calculate_displacement(
            test_args.q0, test_args.v0, mb, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
        )
    with raises(TypeError):
        forced_eqn.calculate_displacement(
            test_args.q0, test_args.v0, 100, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
        )


def test_bad_frequencies(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_eqn.calculate_displacement(
            test_args.q0, test_args.v0, test_args.m, wb, test_args.f, test_args.w, test_args.phi, test_args.t
        )
    with raises(TypeError):
        forced_eqn.calculate_displacement(
            test_args.q0, test_args.v0, test_args.m, 100, test_args.f, test_args.w, test_args.phi, test_args.t
        )
    with raises(errors.UnitsError):
        forced_eqn.calculate_displacement(
            test_args.q0, test_args.v0, test_args.m, test_args.w0, test_args.f, wb, test_args.phi, test_args.t
        )
    with raises(TypeError):
        forced_eqn.calculate_displacement(
            test_args.q0, test_args.v0, test_args.m, test_args.w0, test_args.f, 100, test_args.phi, test_args.t
        )


def test_bad_force(test_args: Args) -> None:
    fb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_eqn.calculate_displacement(
            test_args.q0, test_args.v0, test_args.m, test_args.w0, fb, test_args.w, test_args.phi, test_args.t
        )
    with raises(TypeError):
        forced_eqn.calculate_displacement(
            test_args.q0, test_args.v0, test_args.m, test_args.w0, 100, test_args.w, test_args.phi, test_args.t
        )


def test_bad_phase_lag(test_args: Args) -> None:
    phi_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_eqn.calculate_displacement(
            test_args.q0, test_args.v0, test_args.m, test_args.w0, test_args.f, test_args.w, phi_bad, test_args.t
        )


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_eqn.calculate_displacement(
            test_args.q0, test_args.v0, test_args.m, test_args.w0, test_args.f, test_args.w, test_args.phi, tb
        )
    with raises(TypeError):
        forced_eqn.calculate_displacement(
            test_args.q0, test_args.v0, test_args.m, test_args.w0, test_args.f, test_args.w, test_args.phi, 100
        )
