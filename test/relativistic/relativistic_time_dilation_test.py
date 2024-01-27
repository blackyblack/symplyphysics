from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.relativistic import relativistic_time_dilation


# Using calculations from the paper: https://www.omnicalculator.com/physics/time-dilation


@fixture(name="test_args")
def test_args_fixture():
    own_time = Quantity(10 * units.second)
    velocity = Quantity(200_000_000 * (units.meter / units.second))
    Args = namedtuple("Args", ["own_time", "velocity"])
    return Args(own_time=own_time, velocity=velocity)


def test_basic_time(test_args):
    result = relativistic_time_dilation.calculate_relativistic_time(
        test_args.own_time, test_args.velocity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.time)
    result_own_time = convert_to(result, units.second).evalf(4)
    assert result_own_time == approx(13.423, 0.01)


def test_bad_time(test_args):
    bt = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relativistic_time_dilation.calculate_relativistic_time(bt, test_args.velocity)
    with raises(TypeError):
        relativistic_time_dilation.calculate_relativistic_time(1, test_args.velocity)


def test_bad_velocity(test_args):
    bv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relativistic_time_dilation.calculate_relativistic_time(
            test_args.own_time, bv)
    with raises(TypeError):
        relativistic_time_dilation.calculate_relativistic_time(
            test_args.own_time, 100)
