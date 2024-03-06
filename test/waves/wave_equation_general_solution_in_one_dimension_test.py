from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import wave_equation_general_solution_in_one_dimension as wave_law

# Description
## A wave is traveling in a stretched string with angular frequency w = 5 rad/s and angular
## wavenumber k = 0.5 rad/m, its phase lag is zero. The amplitude of the wave is 10 cm. At
## x = 0.3 m and t = 0.5 s, the displacement of the point of the string is -7.03 cm. Note that
## the cosine form of the solution to the wave equation is assumed.

Args = namedtuple("Args", "a k x w t phi")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(10.0 * units.centimeter)
    k = Quantity(0.5 * units.radian / units.meter)
    x = Quantity(0.3 * units.meter)
    w = Quantity(5.0 * units.radian / units.second)
    t = Quantity(0.5 * units.second)
    phi = 0.0
    return Args(a=a, k=k, x=x, w=w, t=t, phi=phi)


def test_law(test_args: Args) -> None:
    result = wave_law.calculate_displacement(test_args.a, test_args.k, test_args.x, test_args.w, test_args.t, test_args.phi)
    assert_equal(result, -7.03 * units.centimeter)


def test_bad_angular_wavenumber(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wave_law.calculate_displacement(test_args.a, kb, test_args.x, test_args.w, test_args.t, test_args.phi)
    with raises(TypeError):
        wave_law.calculate_displacement(test_args.a, 100, test_args.x, test_args.w, test_args.t, test_args.phi)


def test_bad_position(test_args: Args) -> None:
    xb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wave_law.calculate_displacement(test_args.a, test_args.k, xb, test_args.w, test_args.t, test_args.phi)
    with raises(TypeError):
        wave_law.calculate_displacement(test_args.a, test_args.k, 100, test_args.w, test_args.t, test_args.phi)


def test_bad_angular_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wave_law.calculate_displacement(test_args.a, test_args.k, test_args.x, wb, test_args.t, test_args.phi)
    with raises(TypeError):
        wave_law.calculate_displacement(test_args.a, test_args.k, test_args.x, 100, test_args.t, test_args.phi)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wave_law.calculate_displacement(test_args.a, test_args.k, test_args.x, test_args.w, tb, test_args.phi)
    with raises(TypeError):
        wave_law.calculate_displacement(test_args.a, test_args.k, test_args.x, test_args.w, 100, test_args.phi)


def test_bad_phase_lag(test_args: Args) -> None:
    phi_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wave_law.calculate_displacement(test_args.a, test_args.k, test_args.x, test_args.w, test_args.t, phi_bad)
