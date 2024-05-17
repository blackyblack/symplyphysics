from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.transmission_lines import coefficient_standing_wave_from_voltages as coefficient_law

# Description
## The maximum and minimum voltages in the transmission line are 3 volt and 2 volt, respectively.
## Then the standing wave coefficient is 1.5.

Args = namedtuple("Args", ["maximum_voltage", "minimum_voltage"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    maximum_voltage = Quantity(3 * units.volt)
    minimum_voltage = Quantity(2 * units.volt)

    return Args(maximum_voltage=maximum_voltage, minimum_voltage=minimum_voltage)


def test_basic_coefficient_standing_wave(test_args: Args) -> None:
    result = coefficient_law.calculate_coefficient_standing_wave(test_args.maximum_voltage,
        test_args.minimum_voltage)
    assert_equal(result, 1.5)


def test_bad_voltages(test_args: Args) -> None:
    bad_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_coefficient_standing_wave(bad_voltage, test_args.minimum_voltage)
    with raises(TypeError):
        coefficient_law.calculate_coefficient_standing_wave(100, test_args.minimum_voltage)
    with raises(errors.UnitsError):
        coefficient_law.calculate_coefficient_standing_wave(test_args.maximum_voltage, bad_voltage)
    with raises(TypeError):
        coefficient_law.calculate_coefficient_standing_wave(test_args.maximum_voltage, 100)
    with raises(ValueError):
        coefficient_law.calculate_coefficient_standing_wave(test_args.minimum_voltage,
            test_args.maximum_voltage)
