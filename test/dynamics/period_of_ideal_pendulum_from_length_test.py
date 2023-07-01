from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import period_of_ideal_pendulum_from_length as pendulum_period


@fixture(name="test_args")
def test_args_fixture():
    L = Quantity(1 * units.meter)
    Args = namedtuple("Args", ["L"])
    return Args(L=L)


def test_basic_period(test_args):
    result = pendulum_period.calculate_period(test_args.L)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.time)
    result_period = convert_to(result, units.second).subs(units.second, 1).evalf(2)
    # For a pendulum of 1 meter length, period should be 2 seconds
    assert result_period == approx(2.0, 0.01)


def test_bad_length():
    Lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pendulum_period.calculate_period(Lb)
    with raises(TypeError):
        pendulum_period.calculate_period(100)
