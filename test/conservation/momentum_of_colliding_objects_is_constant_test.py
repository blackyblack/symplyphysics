from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.conservation import momentum_of_colliding_objects_is_constant as conservation_law


@fixture(name="test_args")
def test_args_fixture():
    Ps = Quantity(5 * units.kilogram * units.meter / units.second)
    Args = namedtuple("Args", ["Ps"])
    return Args(Ps=Ps)


def test_basic_conservation(test_args):
    result = conservation_law.calculate_momentum_after(test_args.Ps)
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
        conservation_law.calculate_momentum_after(Pb)
    with raises(AttributeError):
        conservation_law.calculate_momentum_after(100)
