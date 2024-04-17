from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.circuits import wave_resistance_of_lange_coupler as resistance_law

## The wave resistance of the odd mode is 35 ohm, and the wave resistance of the even mode is 25 ohm.
## The number of segments of the Lange coupler is 4. Then the equivalent wave resistance of the coupler is 14.84 ohm.

Args = namedtuple("Args", ["wave_resistance_odd_modes", "wave_resistance_even_modes", "number_segments"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    wave_resistance_odd_modes = Quantity(35 * units.ohm)
    wave_resistance_even_modes = Quantity(25 * units.ohm)
    number_segments = 4
    return Args(wave_resistance_odd_modes=wave_resistance_odd_modes,
        wave_resistance_even_modes=wave_resistance_even_modes,
        number_segments=number_segments)


def test_basic_wave_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_wave_resistance(test_args.wave_resistance_odd_modes,
        test_args.wave_resistance_even_modes, test_args.number_segments)
    assert_equal(result, 14.84 * units.ohm)


def test_bad_resistance(test_args: Args) -> None:
    bad_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance(bad_resistance, test_args.wave_resistance_even_modes,
            test_args.number_segments)
    with raises(TypeError):
        resistance_law.calculate_wave_resistance(100, test_args.wave_resistance_even_modes,
            test_args.number_segments)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance(test_args.wave_resistance_odd_modes, bad_resistance,
            test_args.number_segments)
    with raises(TypeError):
        resistance_law.calculate_wave_resistance(test_args.wave_resistance_odd_modes, 100,
            test_args.number_segments)


def test_bad_number_segments(test_args: Args) -> None:
    bad_number_segments = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance(test_args.wave_resistance_odd_modes,
            test_args.wave_resistance_even_modes, bad_number_segments)
