from collections import namedtuple
from pytest import fixture, raises
from sympy import I
from symplyphysics import (errors, units, Quantity, assert_equal)

from symplyphysics.laws.electricity.circuits.transmission_lines import impedances_t_type_circuit_of_transmission_line as impedances_law

## Characteristic resistance of the transmission line is equal to 50 ohm, line is equal 1 meter,
## constant propogation is equal 6300 [1 / meter], loss factor is equal to 1.7 [1 / meter].
## Then the impedances Z1, Z2, Z3 are equal, respectively: (55.56 - 18.78 * I) ohm,
## (55.56 - 18.78 * I) ohm, (-7.59 + 16.21 * I) ohm.

Args = namedtuple("Args",
    ["characteristic_resistance", "line_length", "constant_propagation"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    characteristic_resistance = Quantity(50 * units.ohm)
    line_length = Quantity(1 * units.meter)
    constant_propagation = Quantity(6300 / units.meter)
    loss_factor = Quantity(1.7 / units.meter)
    actual_propagation_constant = Quantity(loss_factor + I * constant_propagation)
    return Args(characteristic_resistance=characteristic_resistance,
        line_length=line_length,
        constant_propagation=actual_propagation_constant)


def test_basic_transmission_matrix(test_args: Args) -> None:
    result = impedances_law.calculate_impedances(test_args.characteristic_resistance,
        test_args.line_length, test_args.constant_propagation)
    assert_equal(result[0], (55.56 - 18.78 * I) * units.ohm)
    assert_equal(result[1], (55.56 - 18.78 * I) * units.ohm)
    assert_equal(result[2], (-7.59 + 16.21 * I) * units.ohm)


def test_bad_characteristic_resistance(test_args: Args) -> None:
    bad_characteristic_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedances_law.calculate_impedances(bad_characteristic_resistance, test_args.line_length,
            test_args.constant_propagation)
    with raises(TypeError):
        impedances_law.calculate_impedances(100, test_args.line_length,
            test_args.constant_propagation)


def test_bad_line_length(test_args: Args) -> None:
    bad_line_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedances_law.calculate_impedances(test_args.characteristic_resistance, bad_line_length,
            test_args.constant_propagation)
    with raises(TypeError):
        impedances_law.calculate_impedances(test_args.characteristic_resistance, 100,
            test_args.constant_propagation)


def test_bad_constant_propagation(test_args: Args) -> None:
    bad_constant_propagation = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedances_law.calculate_impedances(test_args.characteristic_resistance,
            test_args.line_length, bad_constant_propagation)
    with raises(TypeError):
        impedances_law.calculate_impedances(test_args.characteristic_resistance,
            test_args.line_length, 100)
