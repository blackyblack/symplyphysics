from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import displacement_in_interfering_waves as interference_law

# Description
## Two waves with amplitude u_max = 1 cm, angular wavenumber k = 0.3 rad/m and angular frequency
## w = 50 rad/s are traveling in the same direction and are shifted by pi/2 relative to one another.
## Then at position x = 1 m and time t = 0.5 s the displacement of the resulting wave is 1.33 cm.

Args = namedtuple("Args", "a phi k w x t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(1 * units.centimeter)
    phi = pi / 2
    k = Quantity(0.3 * units.radian / units.meter)
    w = Quantity(50 * units.radian / units.second)
    x = Quantity(1 * units.meter)
    t = Quantity(0.5 * units.second)
    return Args(a=a, phi=phi, k=k, w=w, x=x, t=t)


def test_law(test_args: Args) -> None:
    result = interference_law.calculate_displacement(test_args.a, test_args.phi, test_args.k,
        test_args.w, test_args.x, test_args.t)
    assert_equal(result, 1.33 * units.centimeter, tolerance=3e-3)


def test_bad_phase_shift(test_args: Args) -> None:
    phi_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        interference_law.calculate_displacement(test_args.a, phi_bad, test_args.k, test_args.w,
            test_args.x, test_args.t)


def test_bad_angular_wavenumber(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        interference_law.calculate_displacement(test_args.a, test_args.phi, kb, test_args.w,
            test_args.x, test_args.t)
    with raises(TypeError):
        interference_law.calculate_displacement(test_args.a, test_args.phi, 100, test_args.w,
            test_args.x, test_args.t)


def test_bad_angular_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        interference_law.calculate_displacement(test_args.a, test_args.phi, test_args.k, wb,
            test_args.x, test_args.t)
    with raises(TypeError):
        interference_law.calculate_displacement(test_args.a, test_args.phi, test_args.k, 100,
            test_args.x, test_args.t)


def test_bad_position(test_args: Args) -> None:
    xb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        interference_law.calculate_displacement(test_args.a, test_args.phi, test_args.k,
            test_args.w, xb, test_args.t)
    with raises(TypeError):
        interference_law.calculate_displacement(test_args.a, test_args.phi, test_args.k,
            test_args.w, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        interference_law.calculate_displacement(test_args.a, test_args.phi, test_args.k,
            test_args.w, test_args.x, tb)
    with raises(TypeError):
        interference_law.calculate_displacement(test_args.a, test_args.phi, test_args.k,
            test_args.w, test_args.x, 100)
