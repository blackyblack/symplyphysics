from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.circuits.couplers import relative_operating_band_of_quarter_wave_transformer as bandwidth_law

## Characteristic resistance of the transmission line to which the transformer is connected is equal to 50 ohm.
## Load resistance is equal to 100 ohm. The reflection coefficient is equal to 0.2.

Args = namedtuple("Args",
    ["load_resistance", "characteristic_resistance", "reflection_coefficient"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    load_resistance = Quantity(100 * units.ohm)
    characteristic_resistance = Quantity(50 * units.ohm)
    reflection_coefficient = 0.2
    return Args(load_resistance=load_resistance,
        characteristic_resistance=characteristic_resistance,
        reflection_coefficient=reflection_coefficient)


def test_basic_relative_bandwidth(test_args: Args) -> None:
    result = bandwidth_law.calculate_relative_bandwidth(test_args.load_resistance,
        test_args.characteristic_resistance, test_args.reflection_coefficient)
    assert_equal(result, 0.784)


def test_bad_resistance(test_args: Args) -> None:
    bad_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bandwidth_law.calculate_relative_bandwidth(bad_resistance,
            test_args.characteristic_resistance, test_args.reflection_coefficient)
    with raises(TypeError):
        bandwidth_law.calculate_relative_bandwidth(100, test_args.characteristic_resistance,
            test_args.reflection_coefficient)
    with raises(errors.UnitsError):
        bandwidth_law.calculate_relative_bandwidth(test_args.load_resistance, bad_resistance,
            test_args.reflection_coefficient)
    with raises(TypeError):
        bandwidth_law.calculate_relative_bandwidth(test_args.load_resistance, 100,
            test_args.reflection_coefficient)


def test_bad_reflection_coefficient(test_args: Args) -> None:
    bad_reflection_coefficient = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bandwidth_law.calculate_relative_bandwidth(test_args.load_resistance,
            test_args.characteristic_resistance, bad_reflection_coefficient)
