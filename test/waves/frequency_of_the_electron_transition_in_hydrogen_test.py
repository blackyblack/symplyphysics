from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    units,
    SI,
    Quantity,
    convert_to,
    errors
)
from symplyphysics.laws.waves import frequency_of_the_electron_transition_in_hydrogen as frequency_law

# Description
## It is known that the electron transition frequency is 4.77e14 Hz at m=2 and n=3.
## http://www.phys.utk.edu/labs/modphys/BalmerSeries.pdf

@fixture(name="test_args")
def test_args_fixture():
    number_level_to = 2
    number_level_from = 3

    Args = namedtuple("Args", ["number_level_to", "number_level_from"])
    return Args(number_level_to=number_level_to,
        number_level_from=number_level_from)


def test_basic_transition_frequency(test_args):
    result = frequency_law.calculate_transition_frequency(test_args.number_level_to,
        test_args.number_level_from)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.frequency)
    result_1 = convert_to(result, units.hertz).evalf(5)
    assert result_1 == approx(4.77e14, rel=0.1)


def test_bad_number_level_to(test_args):
    number_level_to = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_transition_frequency(number_level_to, test_args.number_level_from)
    with raises(TypeError):
        frequency_law.calculate_transition_frequency(True, test_args.number_level_from)


def test_bad_number_level_from(test_args):
    number_level_from = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_transition_frequency(test_args.number_level_to, number_level_from)
    with raises(TypeError):
        frequency_law.calculate_transition_frequency(test_args.number_level_to, True)
