from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import resistance_of_resistor_in_wilkinson_divider as resistance_law

# Description
## The standard transmission line resistance is 50 ohms. Then, with the power ratio at the output ports equal to 2,
## the resistance of the resistor will be equal to 125 ohms.

Args = namedtuple("Args", ["transmission_line_resistance", "ratio_of_power"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    transmission_line_resistance = Quantity(50 * units.ohm)
    ratio_of_power = 2

    return Args(transmission_line_resistance=transmission_line_resistance,
        ratio_of_power=ratio_of_power)


def test_basic_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_resistance(test_args.transmission_line_resistance,
        test_args.ratio_of_power)
    assert_equal(result, 125 * units.ohm)


def test_bad_transmission_line_resistance(test_args: Args) -> None:
    transmission_line_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(transmission_line_resistance, test_args.ratio_of_power)
    with raises(TypeError):
        resistance_law.calculate_resistance(100, test_args.ratio_of_power)


def test_bad_ratio_of_power(test_args: Args) -> None:
    ratio_of_power = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(test_args.transmission_line_resistance, ratio_of_power)
