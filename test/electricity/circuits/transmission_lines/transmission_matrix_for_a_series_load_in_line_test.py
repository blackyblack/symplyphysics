from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.circuits.transmission_lines import transmission_matrix_for_a_series_load_in_line as matrix_law

## Load impedance is equal to 50 ohm.
## Then the values of A, B, C, D parameters are equal, respectively: 1, 50 ohm, 0 siemens, 1.

Args = namedtuple("Args", ["load_impedance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    load_impedance = Quantity(50 * units.ohm)
    return Args(load_impedance=load_impedance)


def test_basic_transmission_matrix(test_args: Args) -> None:
    result = matrix_law.calculate_transmission_matrix(test_args.load_impedance)
    assert_equal(result[0][0], 1)
    assert_equal(result[0][1], 50 * units.ohm)
    assert_equal(result[1][0], 0 * units.siemens)
    assert_equal(result[1][1], 1)


def test_bad_load_impedance() -> None:
    bad_load_impedance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(bad_load_impedance)
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(100)
