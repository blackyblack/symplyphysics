from collections import namedtuple
from pytest import fixture, raises
from sympy import I
from symplyphysics import (errors, units, Quantity, assert_equal)

from symplyphysics.laws.electricity.circuits.transmission_lines import input_impedance_lossy_transmission_line as impedance_law

## Characteristic resistance of the transmission line is equal to 50 ohm, line is equal 1 meter,
## constant propogation is equal 6300 [1 / meter], loss factor is equal to 1.7 [1 / meter].
## Load resistance is equal to 100 ohm.
## Then the input impedance is equal to 49.33 - 0.879 * I ohm.

Args = namedtuple("Args", [
    "characteristic_resistance", "load_resistance", "constant_propagation", "line_length",
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    characteristic_resistance = Quantity(50 * units.ohm)
    load_resistance = Quantity(100 * units.ohm)
    constant_propagation = Quantity(6300 / units.meter)
    line_length = Quantity(1 * units.meter)
    loss_factor = Quantity(1.7 / units.meter)
    actual_propagation_constant = Quantity(loss_factor + I * constant_propagation)
    return Args(
        characteristic_resistance=characteristic_resistance,
        load_resistance=load_resistance,
        constant_propagation=actual_propagation_constant,
        line_length=line_length,
    )


def test_basic_input_impedance(test_args: Args) -> None:
    result = impedance_law.calculate_input_impedance(test_args.characteristic_resistance,
        test_args.load_resistance, test_args.constant_propagation, test_args.line_length)
    assert_equal(result, (49.33 - 0.879 * I) * units.ohm)


def test_bad_characteristic_resistance(test_args: Args) -> None:
    bad_characteristic_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedance_law.calculate_input_impedance(bad_characteristic_resistance,
            test_args.load_resistance, test_args.constant_propagation, test_args.line_length)
    with raises(TypeError):
        impedance_law.calculate_input_impedance(100, test_args.load_resistance,
            test_args.constant_propagation, test_args.line_length)


def test_bad_load_resistance(test_args: Args) -> None:
    bad_load_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedance_law.calculate_input_impedance(test_args.characteristic_resistance,
            bad_load_resistance, test_args.constant_propagation, test_args.line_length)
    with raises(TypeError):
        impedance_law.calculate_input_impedance(test_args.characteristic_resistance, 100,
            test_args.constant_propagation, test_args.line_length)


def test_bad_constant_propagation(test_args: Args) -> None:
    constant_propagation = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedance_law.calculate_input_impedance(test_args.characteristic_resistance,
            test_args.load_resistance, constant_propagation, test_args.line_length)
    with raises(TypeError):
        impedance_law.calculate_input_impedance(test_args.characteristic_resistance,
            test_args.load_resistance, 100, test_args.line_length)


def test_bad_line_length(test_args: Args) -> None:
    line_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedance_law.calculate_input_impedance(test_args.characteristic_resistance,
            test_args.load_resistance, test_args.constant_propagation, line_length)
    with raises(TypeError):
        impedance_law.calculate_input_impedance(test_args.characteristic_resistance,
            test_args.load_resistance, test_args.constant_propagation, 100)
