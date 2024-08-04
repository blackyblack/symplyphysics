from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.circuits.couplers import wave_resistance_even_mode_of_lange_coupler as resistance_law

## The wave resistance of the odd mode is 80 ohm. The coupling coefficient between the segments of the
## coupler is 0.3. The number of segments of the Lange coupler is 4. Then the wave resistance of the even mode is 121.05 ohm.

Args = namedtuple("Args", ["coupling_factor", "wave_resistance_odd_modes", "number_segments"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    coupling_factor = 0.3
    wave_resistance_odd_modes = Quantity(80 * units.ohm)
    number_segments = 4
    return Args(coupling_factor=coupling_factor,
        wave_resistance_odd_modes=wave_resistance_odd_modes,
        number_segments=number_segments)


def test_basic_wave_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_wave_resistance_even_modes(
        test_args.coupling_factor, test_args.wave_resistance_odd_modes, test_args.number_segments)
    assert_equal(result, 121.05 * units.ohm)


def test_bad_coupling_factor(test_args: Args) -> None:
    bad_coupling_factor = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance_even_modes(bad_coupling_factor,
            test_args.wave_resistance_odd_modes, test_args.number_segments)


def test_bad_wave_resistance_odd_modes(test_args: Args) -> None:
    bad_wave_resistance_odd_modes = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance_even_modes(test_args.coupling_factor,
            bad_wave_resistance_odd_modes, test_args.number_segments)
    with raises(TypeError):
        resistance_law.calculate_wave_resistance_even_modes(test_args.coupling_factor, 100,
            test_args.number_segments)


def test_bad_number_segments(test_args: Args) -> None:
    bad_number_segments = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance_even_modes(test_args.coupling_factor,
            test_args.wave_resistance_odd_modes, bad_number_segments)
