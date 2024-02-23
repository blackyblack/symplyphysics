from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import (
    wave_speed_from_angular_frequency_and_wavenumber as wave_speed_law,
)

# Description
## A wave has angular frequency w = 9 Hz and angular wavenumber k = 0.3 1/m.
## Its speed is 30 m/s.

Args = namedtuple("Args", "w k")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w = Quantity(9.0 * units.hertz)
    k = Quantity(0.3 / units.meter)
    return Args(w=w, k=k)


def test_law(test_args: Args) -> None:
    result = wave_speed_law.calculate_wave_speed(test_args.w, test_args.k)
    assert_equal(result, 30.0 * units.meter / units.second)


def test_bad_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wave_speed_law.calculate_wave_speed(wb, test_args.k)
    with raises(TypeError):
        wave_speed_law.calculate_wave_speed(100, test_args.k)


def test_bad_wavenumber(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wave_speed_law.calculate_wave_speed(test_args.w, kb)
    with raises(TypeError):
        wave_speed_law.calculate_wave_speed(test_args.w, 100)
