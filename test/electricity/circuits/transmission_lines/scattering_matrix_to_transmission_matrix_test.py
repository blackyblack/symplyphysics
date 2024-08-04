from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.circuits.transmission_lines import scattering_matrix_to_transmission_matrix as matrix_law

## S-parameters S11, S12, S21, S22 are equal to 0.2, 0.5, 0.5, 0.2, respectively.
## Characteristic resistance of the transmission line is equal to 50 ohm.
## Then the values of A, B, C, D parameters are equal, respectively: 1.21, 59.5 ohm, 7.8 millisiemens, 1.21.

Args = namedtuple("Args", ["characteristic_resistance", "parameters"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    characteristic_resistance = Quantity(50 * units.ohm)
    parameters = ((0.2, 0.5), (0.5, 0.2))
    return Args(characteristic_resistance=characteristic_resistance, parameters=parameters)


def test_basic_waves(test_args: Args) -> None:
    result = matrix_law.calculate_transmission_matrix(test_args.characteristic_resistance,
        test_args.parameters)
    assert_equal(result[0][0], 1.21)
    assert_equal(result[0][1], 59.5 * units.ohm)
    assert_equal(result[1][0], 7.8 * prefixes.milli * units.siemens)
    assert_equal(result[1][1], 1.21)


def test_bad_characteristic_resistance(test_args: Args) -> None:
    bad_characteristic_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(bad_characteristic_resistance,
            test_args.parameters)
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(100, test_args.parameters)
