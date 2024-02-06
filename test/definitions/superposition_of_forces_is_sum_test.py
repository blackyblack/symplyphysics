from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import superposition_of_forces_is_sum as forces_law


@fixture(name="test_args")
def test_args_fixture():
    F1 = Quantity(10 * units.newton)
    F2 = Quantity(20 * units.newton)
    Args = namedtuple("Args", ["F1", "F2"])
    return Args(F1=F1, F2=F2)


def test_basic_superposition(test_args):
    result = forces_law.calculate_resultant_force([test_args.F1, test_args.F2])
    assert_equal(result, 30 * units.newton)


def test_three_forces_array(test_args):
    F3 = Quantity(-5 * units.newton)
    result = forces_law.calculate_resultant_force([test_args.F1, test_args.F2, F3])
    assert_equal(result, 25 * units.newton)


def test_bad_force(test_args):
    Fb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([Fb, test_args.F2])
    with raises(TypeError):
        forces_law.calculate_resultant_force([100, test_args.F2])
    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([test_args.F1, Fb])
    with raises(TypeError):
        forces_law.calculate_resultant_force([test_args.F1, 100])
    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([Fb, Fb])
    with raises(TypeError):
        forces_law.calculate_resultant_force([100, 100])
    with raises(TypeError):
        forces_law.calculate_resultant_force(test_args.F1)
