from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.conservation import momentum_of_colliding_objects_is_constant as conservation_law


@fixture(name="test_args")
def test_args_fixture():
    Ps = Quantity(5 * units.kilogram * units.meter / units.second)
    Args = namedtuple("Args", ["Ps"])
    return Args(Ps=Ps)


def test_basic_conservation(test_args):
    result = conservation_law.calculate_momentum_after(test_args.Ps)
    assert_equal(result, 5 * units.kilogram * units.meter / units.second)


def test_bad_momentum():
    Pb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        conservation_law.calculate_momentum_after(Pb)
    with raises(TypeError):
        conservation_law.calculate_momentum_after(100)
