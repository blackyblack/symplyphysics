from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.transmission_lines import transmission_matrix_of_t_type_circuit as matrix_law

## First impedance is equal to 50 ohm, second impedance is equal to 50 ohm,
## third impedance is equal to 100 ohm.
## Then the values of A, B, C, D parameters are equal, respectively: 1.5, 125 ohm, 10 millisiemens, 1.5.

Args = namedtuple("Args", ["first_impedance", "second_impedance", "third_impedance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    first_impedance = Quantity(50 * units.ohm)
    second_impedance = Quantity(50 * units.ohm)
    third_impedance = Quantity(100 * units.ohm)
    return Args(
        first_impedance=first_impedance,
        second_impedance=second_impedance,
        third_impedance=third_impedance,
    )


def test_basic_transmission_matrix(test_args: Args) -> None:
    result = matrix_law.calculate_transmission_matrix(
        (test_args.first_impedance, test_args.second_impedance, test_args.third_impedance))
    assert_equal(result[0][0], 1.5)
    assert_equal(result[0][1], 125 * units.ohm)
    assert_equal(result[1][0], 10 * prefixes.milli * units.siemens)
    assert_equal(result[1][1], 1.5)


def test_bad_impedances(test_args: Args) -> None:
    bad_impedance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(
            (bad_impedance, test_args.second_impedance, test_args.third_impedance))
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(
            (100, test_args.second_impedance, test_args.third_impedance))
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(
            (test_args.first_impedance, bad_impedance, test_args.third_impedance))
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(
            (test_args.first_impedance, 100, test_args.third_impedance))
    with raises(errors.UnitsError):
        matrix_law.calculate_transmission_matrix(
            (test_args.first_impedance, test_args.second_impedance, bad_impedance))
    with raises(TypeError):
        matrix_law.calculate_transmission_matrix(
            (test_args.first_impedance, test_args.second_impedance, 100))
