from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import momentum_after_collision_equals_to_momentum_before as momentum_law


@fixture
def test_args():
    P_before = Quantity(5 * units.kilogram * units.meter / units.second)
    Args = namedtuple("Args", ["P_before"])
    return Args(P_before=P_before)


def test_basic_conservation(test_args):
    result = momentum_law.calculate_momentum_after(test_args.P_before)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.momentum)
    result_ = convert_to(result, units.kilogram * units.meter / units.second).subs({
        units.kilogram: 1,
        units.meter: 1,
        units.second: 1
    }).evalf(2)
    assert result_ == approx(5.0, 0.01)


def test_bad_momentum():
    Pb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        momentum_law.calculate_momentum_after(Pb)
    with raises(TypeError):
        momentum_law.calculate_momentum_after(100)
