from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import resonant_oscillations as resonance_law

# Description
## An external driving force is acting on an oscillator and is in resonance with it.
## The amplitude of the force is 1 N and its phase lag is zero. The oscillator's mass
## is 0.5 kg and its natural angular frequency is 10 rad/s. At t = 10 s, the displacement
## of the oscillator due to the driving force amounts to -50.6 cm.

Args = namedtuple("Args", "m w0 f phi t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(0.5 * units.kilogram)
    w0 = Quantity(10.0 * units.radian / units.second)
    f = Quantity(1.0 * units.newton)
    phi = 0.0
    t = Quantity(10.0 * units.second)
    return Args(m=m, w0=w0, f=f, phi=phi, t=t)


def test_law(test_args: Args) -> None:
    result = resonance_law.calculate_resonant_displacement(test_args.m, test_args.w0, test_args.f, test_args.phi, test_args.t)
    assert_equal(result, -50.6 * units.centimeter)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        resonance_law.calculate_resonant_displacement(mb, test_args.w0, test_args.f, test_args.phi, test_args.t)
    with raises(TypeError):
        resonance_law.calculate_resonant_displacement(100, test_args.w0, test_args.f, test_args.phi, test_args.t)


def test_bad_natural_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        resonance_law.calculate_resonant_displacement(test_args.m, wb, test_args.f, test_args.phi, test_args.t)
    with raises(TypeError):
        resonance_law.calculate_resonant_displacement(test_args.m, 100, test_args.f, test_args.phi, test_args.t)


def test_bad_force_amplitude(test_args: Args) -> None:
    fb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        resonance_law.calculate_resonant_displacement(test_args.m, test_args.w0, fb, test_args.phi, test_args.t)
    with raises(TypeError):
        resonance_law.calculate_resonant_displacement(test_args.m, test_args.w0, 100, test_args.phi, test_args.t)


def test_bad_phase_lag(test_args: Args) -> None:
    phi_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        resonance_law.calculate_resonant_displacement(test_args.m, test_args.w0, test_args.f, phi_bad, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        resonance_law.calculate_resonant_displacement(test_args.m, test_args.w0, test_args.f, test_args.phi, tb)
    with raises(TypeError):
        resonance_law.calculate_resonant_displacement(test_args.m, test_args.w0, test_args.f, test_args.phi, 100)
