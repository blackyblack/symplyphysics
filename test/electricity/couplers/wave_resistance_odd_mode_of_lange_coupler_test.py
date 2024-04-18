from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.couplers import wave_resistance_odd_mode_of_lange_coupler as resistance_law

## The standard characteristic resistance of the transmission line is 50 ohm. The coupling coefficient between the
## segments of the coupler is 0.3. The number of segments of the Lange coupler is 4. Then the wave resistance of the odd mode is 80.87 ohm.

Args = namedtuple("Args", ["coupling_factor", "characteristic_resistance", "number_segments"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    coupling_factor = 0.3
    characteristic_resistance = Quantity(50 * units.ohm)
    number_segments = 4
    return Args(coupling_factor=coupling_factor,
        characteristic_resistance=characteristic_resistance,
        number_segments=number_segments)


def test_basic_wave_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_wave_resistance_odd_modes(test_args.coupling_factor,
        test_args.characteristic_resistance, test_args.number_segments)
    assert_equal(result, 80.87 * units.ohm)


def test_bad_coupling_factor(test_args: Args) -> None:
    bad_coupling_factor = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance_odd_modes(bad_coupling_factor, test_args.characteristic_resistance,
            test_args.number_segments)


def test_bad_characteristic_resistance(test_args: Args) -> None:
    bad_characteristic_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance_odd_modes(test_args.coupling_factor, bad_characteristic_resistance,
            test_args.number_segments)
    with raises(TypeError):
        resistance_law.calculate_wave_resistance_odd_modes(test_args.coupling_factor, 100,
            test_args.number_segments)


def test_bad_number_segments(test_args: Args) -> None:
    bad_number_segments = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance_odd_modes(test_args.coupling_factor,
            test_args.characteristic_resistance, bad_number_segments)
