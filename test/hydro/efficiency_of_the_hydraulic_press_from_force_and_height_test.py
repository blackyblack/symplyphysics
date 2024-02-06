from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.laws.hydro import efficiency_of_the_hydraulic_press_from_force_and_height as efficiency


@fixture(name="test_args")
def test_args_fixture():
    useful_force = Quantity(9500 * units.newton)
    useful_height = Quantity(0.01 * units.meter)
    expended_force = Quantity(500 * units.newton)
    expended_height = Quantity(0.2 * units.meter)
    Args = namedtuple("Args",
        ["useful_force", "useful_height", "expended_force", "expended_height"])
    return Args(useful_force=useful_force,
        useful_height=useful_height,
        expended_force=expended_force,
        expended_height=expended_height)


def test_basic_efficiency(test_args):
    result = efficiency.calculate_efficiency(test_args.useful_force, test_args.useful_height,
        test_args.expended_force, test_args.expended_height)
    assert_equal(result, 0.95)


def test_bad_force(test_args):
    bad_force = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        efficiency.calculate_efficiency(bad_force, test_args.useful_height,
            test_args.expended_force, test_args.expended_height)
    with raises(errors.UnitsError):
        efficiency.calculate_efficiency(test_args.useful_force, test_args.useful_height, bad_force,
            test_args.expended_height)
    with raises(TypeError):
        efficiency.calculate_efficiency(100, test_args.useful_height, test_args.expended_force,
            test_args.expended_height)
    with raises(TypeError):
        efficiency.calculate_efficiency(test_args.useful_force, test_args.useful_height, 100,
            test_args.expended_height)


def test_bad_height(test_args):
    bad_height = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        efficiency.calculate_efficiency(test_args.useful_force, bad_height,
            test_args.expended_force, test_args.expended_height)
    with raises(errors.UnitsError):
        efficiency.calculate_efficiency(test_args.useful_force, test_args.useful_height,
            test_args.expended_force, bad_height)
    with raises(TypeError):
        efficiency.calculate_efficiency(test_args.useful_force, 100, test_args.expended_force,
            test_args.expended_height)
    with raises(TypeError):
        efficiency.calculate_efficiency(test_args.useful_force, test_args.useful_height,
            test_args.expended_force, 100)
