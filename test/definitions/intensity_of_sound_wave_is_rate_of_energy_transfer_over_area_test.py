from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.definitions import (
    intensity_of_sound_wave_is_rate_of_energy_transfer_over_area as intensity_def,
)

# Description
## A sound wave of power P = 1 pW is hitting a surface of area A = 0.5 m**2. Then its
## intensity is 2 pW/m**2.

Args = namedtuple("Args", "p a")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(1 * prefixes.pico * units.watt)
    a = Quantity(0.5 * units.meter**2)
    return Args(p=p, a=a)


def test_law(test_args: Args) -> None:
    result = intensity_def.calculate_intensity(test_args.p, test_args.a)
    assert_equal(result, 2 * prefixes.pico * units.watt / units.meter**2)


def test_bad_power(test_args: Args) -> None:
    pb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_def.calculate_intensity(pb, test_args.a)
    with raises(TypeError):
        intensity_def.calculate_intensity(100, test_args.a)


def test_bad_area(test_args: Args) -> None:
    ab = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_def.calculate_intensity(test_args.p, ab)
    with raises(TypeError):
        intensity_def.calculate_intensity(test_args.p, 100)
