from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal)

from symplyphysics.laws.electricity.couplers import impedances_for_wilkinson_microstrip_divider as impedances_law

## The resistance of the transmission line to which the divider is connected is 50 ohms, the power ratio at the output ports equal to 2.
## Then the values of Z1, Z2, Z3, Z4 impedances are equal, respectively: 158.11 ohm, 39.53 ohm, 70.71 ohm, 35.36 ohm.

Args = namedtuple("Args", ["characteristic_resistance", "ratio_of_power"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    characteristic_resistance = Quantity(50 * units.ohm)
    ratio_of_power = 2

    return Args(characteristic_resistance=characteristic_resistance, ratio_of_power=ratio_of_power)


def test_basic_impedances(test_args: Args) -> None:
    result = impedances_law.calculate_impedances(test_args.characteristic_resistance,
        test_args.ratio_of_power)
    assert_equal(result[0], 158.11 * units.ohm)
    assert_equal(result[1], 39.53 * units.ohm)
    assert_equal(result[2], 70.71 * units.ohm)
    assert_equal(result[3], 35.36 * units.ohm)


def test_bad_characteristic_resistance(test_args: Args) -> None:
    bad_characteristic_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedances_law.calculate_impedances(bad_characteristic_resistance, test_args.ratio_of_power)
    with raises(TypeError):
        impedances_law.calculate_impedances(100, test_args.ratio_of_power)


def test_bad_ratio_of_power(test_args: Args) -> None:
    bad_ratio_of_power = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedances_law.calculate_impedances(test_args.characteristic_resistance, bad_ratio_of_power)
