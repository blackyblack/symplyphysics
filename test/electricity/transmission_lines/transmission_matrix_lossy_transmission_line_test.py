from collections import namedtuple
from pytest import fixture, raises
from sympy import I
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.transmission_lines import transmission_matrix_lossy_transmission_line as matrix_law

## Characteristic resistance of the transmission line is equal to 50 ohm, line is equal 1 meter,
## constant propogation is equal 6300 [1 / meter], loss factor is equal to 1.7 [1 / meter].
## Then the values of A, B, C, D parameters are equal, respectively: -1.265 - 2.365 * I,
## -59.2 - 126.5 * I ohm, -23.7 - 50.6 * I millisiemens, -1.265 - 2.365 * I.

Args = namedtuple("Args",
    ["characteristic_resistance", "line_length", "constant_propagation", "loss_factor"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    characteristic_resistance = Quantity(50 * units.ohm)
    line_length = Quantity(1 * units.meter)
    constant_propagation = Quantity(6300 / units.meter)
    loss_factor = Quantity(1.7 / units.meter)
    return Args(characteristic_resistance=characteristic_resistance,
        line_length=line_length,
        constant_propagation=constant_propagation,
        loss_factor=loss_factor)


def test_basic_transmission_matrix(test_args: Args) -> None:
    result = matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance,
        test_args.line_length, test_args.constant_propagation, test_args.loss_factor)
    assert_equal(result[0][0], -1.265 - 2.365 * I)
    assert_equal(result[0][1], (-59.2 - 126.5 * I) * units.ohm)
    assert_equal(result[1][0], (-23.7 - 50.6 * I) * prefixes.milli * units.siemens)
    assert_equal(result[1][1], -1.265 - 2.365 * I)


def test_bad_characteristic_resistance(test_args: Args) -> None:
    bad_characteristic_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(bad_characteristic_resistance,
            test_args.line_length, test_args.constant_propagation, test_args.loss_factor)
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(100, test_args.line_length,
            test_args.constant_propagation, test_args.loss_factor)


def test_bad_line_length(test_args: Args) -> None:
    bad_line_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance,
            bad_line_length, test_args.constant_propagation, test_args.loss_factor)
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance, 100,
            test_args.constant_propagation, test_args.loss_factor)


def test_bad_constant_propagation(test_args: Args) -> None:
    bad_constant_propagation = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance,
            test_args.line_length, bad_constant_propagation, test_args.loss_factor)
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance,
            test_args.line_length, 100, test_args.loss_factor)


def test_bad_loss_factor(test_args: Args) -> None:
    bad_loss_factor = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance,
            test_args.line_length, test_args.constant_propagation, bad_loss_factor)
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance,
            test_args.line_length, test_args.constant_propagation, 100)
