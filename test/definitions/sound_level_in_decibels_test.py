from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import sound_level_in_decibels as sound_level_def

# Description
## The sound level of a sound wave of intensity I = 1e-3 W/m**2 is 90 dB.

Args = namedtuple("Args", "i")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    i = Quantity(1e-3 * units.watt / units.meter**2)
    return Args(i=i)


def test_law(test_args: Args) -> None:
    result = sound_level_def.calculate_sound_level(test_args.i)
    assert_equal(result, 90)


def test_bad_intensity() -> None:
    ib = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        sound_level_def.calculate_sound_level(ib)
    with raises(TypeError):
        sound_level_def.calculate_sound_level(100)
