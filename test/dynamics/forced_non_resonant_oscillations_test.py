from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import forced_non_resonant_oscillations as forced_law

# Description
## An oscillating external driving force of magnitude f = 1 N and angular frequency w = 5 N is acting on
## an oscillator with a natural angular frequency w0 = 3 N. At t = 1 s, the displacement of the oscillator
## caused by the driving force amounts to q = -8.87 mm. The mass of the oscillator is 2 kg.

Args = namedtuple("Args", "m w0 f w phi t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(2.0 * units.kilogram)
    w0 = Quantity(3.0 * units.radian / units.second)
    f = Quantity(1.0 * units.newton)
    w = Quantity(5.0 * units.radian / units.second)
    phi = 0.0
    t = Quantity(1.0 * units.second)
    return Args(m=m, w0=w0, f=f, w=w, phi=phi, t=t)


def test_law(test_args: Args) -> None:
    result = forced_law.calculate_driven_displacement(
        test_args.m, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
    )
    assert_equal(result, -8.87 * units.millimeter)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_law.calculate_driven_displacement(
            mb, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
        )
    with raises(TypeError):
        forced_law.calculate_driven_displacement(
            100, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
        )


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_law.calculate_driven_displacement(
            mb, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
        )
    with raises(TypeError):
        forced_law.calculate_driven_displacement(
            100, test_args.w0, test_args.f, test_args.w, test_args.phi, test_args.t
        )


def test_bad_angular_frequencies(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_law.calculate_driven_displacement(
            test_args.m, wb, test_args.f, test_args.w, test_args.phi, test_args.t
        )
    with raises(TypeError):
        forced_law.calculate_driven_displacement(
            test_args.m, 100, test_args.f, test_args.w, test_args.phi, test_args.t
        )
    with raises(errors.UnitsError):
        forced_law.calculate_driven_displacement(
            test_args.m, test_args.w0, test_args.f, wb, test_args.phi, test_args.t
        )
    with raises(TypeError):
        forced_law.calculate_driven_displacement(
            test_args.m, test_args.w0, test_args.f, 100, test_args.phi, test_args.t
        )


def test_bad_force_amplitude(test_args: Args) -> None:
    fb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_law.calculate_driven_displacement(
            test_args.m, test_args.w0, fb, test_args.w, test_args.phi, test_args.t
        )
    with raises(TypeError):
        forced_law.calculate_driven_displacement(
            test_args.m, test_args.w0, 100, test_args.w, test_args.phi, test_args.t
        )


def test_bad_phase_lag(test_args: Args) -> None:
    phi_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_law.calculate_driven_displacement(
            test_args.m, test_args.w0, test_args.f, test_args.w, phi_bad, test_args.t
        )


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        forced_law.calculate_driven_displacement(
            test_args.m, test_args.w0, test_args.f, test_args.w, test_args.phi, tb
        )
    with raises(TypeError):
        forced_law.calculate_driven_displacement(
            test_args.m, test_args.w0, test_args.f, test_args.w, test_args.phi, 100
        )
