from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import displacement_in_standing_wave as standing_wave_law

# Description
## A standing wave is occurring on a stretched string with fixed ends. The amplitude of the traveling
## waves comprising it is 10.0 cm. The angular wavenumber of the waves is 5.0 rad/m and the angular frequency
## is 1.0 rad/s. At x = 9.0 cm and t = 2.0 s, the displacement of the string is -3.6 cm.

Args = namedtuple("Args", "a k x w t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(10.0 * units.centimeter)
    k = Quantity(5.0 * units.radian / units.meter)
    x = Quantity(9.0 * units.centimeter)
    w = Quantity(1.0 * units.radian / units.second)
    t = Quantity(2.0 * units.second)
    return Args(a=a, k=k, x=x, w=w, t=t)


def test_law(test_args: Args) -> None:
    result = standing_wave_law.calculate_displacement(test_args.a, test_args.k, test_args.x,
        test_args.w, test_args.t)
    assert_equal(result, -3.6 * units.centimeter, relative_tolerance=6e-3)


def test_bad_angular_wavenumber(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        standing_wave_law.calculate_displacement(test_args.a, kb, test_args.x, test_args.w,
            test_args.t)
    with raises(TypeError):
        standing_wave_law.calculate_displacement(test_args.a, 100, test_args.x, test_args.w,
            test_args.t)


def test_bad_position(test_args: Args) -> None:
    xb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        standing_wave_law.calculate_displacement(test_args.a, test_args.k, xb, test_args.w,
            test_args.t)
    with raises(TypeError):
        standing_wave_law.calculate_displacement(test_args.a, test_args.k, 100, test_args.w,
            test_args.t)


def test_bad_angular_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        standing_wave_law.calculate_displacement(test_args.a, test_args.k, test_args.x, wb,
            test_args.t)
    with raises(TypeError):
        standing_wave_law.calculate_displacement(test_args.a, test_args.k, test_args.x, 100,
            test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        standing_wave_law.calculate_displacement(test_args.a, test_args.k, test_args.x, test_args.w,
            tb)
    with raises(TypeError):
        standing_wave_law.calculate_displacement(test_args.a, test_args.k, test_args.x, test_args.w,
            100)
