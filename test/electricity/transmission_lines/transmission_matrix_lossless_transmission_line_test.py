from collections import namedtuple
from pytest import fixture, raises
from sympy import I
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.transmission_lines import transmission_matrix_lossless_transmission_line as matrix_law

## Characteristic resistance of the transmission line is equal to 50 ohm, line is equal 1 meter,
## constant propogation is equal 6300 [1 / meter].
## Then the values of A, B, C, D parameters are equal, respectively: -0.4476, -44.7 * I ohm, -17.88 * I millisiemens, -0.4476.

Args = namedtuple("Args", ["characteristic_resistance", "line_length", "constant_propagation"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    characteristic_resistance = Quantity(50 * units.ohm)
    line_length = Quantity(1 * units.meter)
    constant_propagation = Quantity(6300 / units.meter)
    return Args(characteristic_resistance=characteristic_resistance,
        line_length=line_length,
        constant_propagation=constant_propagation,
        )


def test_basic_transmission_matrix(test_args: Args) -> None:
    result = matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance, test_args.line_length, test_args.constant_propagation)
    assert_equal(result[0][0], -0.4476)
    assert_equal(result[0][1], -44.7 * I * units.ohm)
    assert_equal(result[1][0], -17.88 * I * prefixes.milli * units.siemens)
    assert_equal(result[1][1], -0.4476)


def test_bad_characteristic_resistance(test_args: Args) -> None:
    bad_characteristic_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(bad_characteristic_resistance, test_args.line_length, test_args.constant_propagation)
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(100, test_args.line_length, test_args.constant_propagation)


def test_bad_line_length(test_args: Args) -> None:
    bad_line_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance, bad_line_length, test_args.constant_propagation)
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance, 100, test_args.constant_propagation)


def test_bad_constant_propagation(test_args: Args) -> None:
    bad_constant_propagation = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance, test_args.line_length, bad_constant_propagation)
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance, test_args.line_length, 100)
