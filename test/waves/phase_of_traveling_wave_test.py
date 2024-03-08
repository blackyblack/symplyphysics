from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import phase_of_traveling_wave as phase_law

# Description
## A traveling wave has angular wavenumber k = 3 rad/m and angular frequency w = 0.8 Hz.
## Then at position x = 0.1 m and time t = 10 s its phase is -7.7 rad. The wave is moving
## in the positive direction of the x-axis.

Args = namedtuple("Args", "k x w t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    k = Quantity(3.0 * units.radian / units.meter)
    x = Quantity(0.1 * units.meter)
    w = Quantity(0.8 * units.hertz)
    t = Quantity(10.0 * units.second)
    return Args(k=k, x=x, w=w, t=t)


def test_law(test_args: Args) -> None:
    result = phase_law.calculate_wave_phase(test_args.k, test_args.x, test_args.w, test_args.t)
    assert_equal(result, -7.7 * units.radian)


def test_bad_wavenumber(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        phase_law.calculate_wave_phase(kb, test_args.x, test_args.w, test_args.t)
    with raises(TypeError):
        phase_law.calculate_wave_phase(100, test_args.x, test_args.w, test_args.t)


def test_bad_position(test_args: Args) -> None:
    xb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        phase_law.calculate_wave_phase(test_args.k, xb, test_args.w, test_args.t)
    with raises(TypeError):
        phase_law.calculate_wave_phase(test_args.k, 100, test_args.w, test_args.t)


def test_bad_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        phase_law.calculate_wave_phase(test_args.k, test_args.x, wb, test_args.t)
    with raises(TypeError):
        phase_law.calculate_wave_phase(test_args.k, test_args.x, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        phase_law.calculate_wave_phase(test_args.k, test_args.x, test_args.w, tb)
    with raises(TypeError):
        phase_law.calculate_wave_phase(test_args.k, test_args.x, test_args.w, 100)
