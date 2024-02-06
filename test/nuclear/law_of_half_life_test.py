from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_approx, units, Quantity, errors)
from symplyphysics.laws.nuclear import law_of_half_life as number_of_cores_law

# Description
## It is known that for oxygen, the half-life is 124 seconds.
## https://ru.wikipedia.org/wiki/Период_полураспада


@fixture(name="test_args")
def test_args_fixture():
    number_of_cores_initial = 2e30
    half_life = Quantity(90 * units.second)
    decay_time = Quantity(90 * units.second)

    Args = namedtuple("Args", ["number_of_cores_initial", "half_life", "decay_time"])
    return Args(number_of_cores_initial=number_of_cores_initial,
        half_life=half_life,
        decay_time=decay_time)


def test_basic_number_of_cores(test_args):
    result = number_of_cores_law.calculate_number_of_cores(test_args.number_of_cores_initial,
        test_args.half_life, test_args.decay_time)
    assert_approx(result, 1e30)


def test_bad_number_of_cores_initial(test_args):
    number_of_cores_initial = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        number_of_cores_law.calculate_number_of_cores(number_of_cores_initial, test_args.half_life,
            test_args.decay_time)
    with raises(ValueError):
        number_of_cores_law.calculate_number_of_cores(-1, test_args.half_life, test_args.decay_time)


def test_bad_half_life(test_args):
    half_life = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        number_of_cores_law.calculate_number_of_cores(test_args.number_of_cores_initial, half_life,
            test_args.decay_time)
    with raises(TypeError):
        number_of_cores_law.calculate_number_of_cores(test_args.number_of_cores_initial, 100,
            test_args.decay_time)


def test_bad_decay_time(test_args):
    decay_time = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        number_of_cores_law.calculate_number_of_cores(test_args.number_of_cores_initial,
            test_args.half_life, decay_time)
    with raises(TypeError):
        number_of_cores_law.calculate_number_of_cores(test_args.number_of_cores_initial,
            test_args.half_life, 100)
