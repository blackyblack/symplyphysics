from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.waves import frequency_of_electron_transition_in_hydrogen as frequency_law

# Description
## It is known that the electron transition frequency is 4.56e14 Hz at m=2 and n=3.
## http://www.phys.utk.edu/labs/modphys/BalmerSeries.pdf

Args = namedtuple("Args", ["number_level_to", "number_level_from"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    number_level_to = 2
    number_level_from = 3
    return Args(number_level_to=number_level_to, number_level_from=number_level_from)


def test_basic_transition_frequency(test_args: Args) -> None:
    result = frequency_law.calculate_transition_frequency(test_args.number_level_to,
        test_args.number_level_from)
    assert_equal(result, 4.568e14 * units.hertz)


def test_bad_number_level_to(test_args: Args) -> None:
    number_level_to = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_transition_frequency(number_level_to, test_args.number_level_from)
    with raises(TypeError):
        frequency_law.calculate_transition_frequency(True, test_args.number_level_from)


def test_bad_number_level_from(test_args: Args) -> None:
    number_level_from = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_transition_frequency(test_args.number_level_to, number_level_from)
    with raises(TypeError):
        frequency_law.calculate_transition_frequency(test_args.number_level_to, True)
