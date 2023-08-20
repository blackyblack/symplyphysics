from collections import namedtuple
from pytest import approx, fixture
from symplyphysics import (
    units,
    convert_to,
    Quantity,
    SI,
    QuantityVector,
)
from symplyphysics.laws.dynamics.vector import acceleration_from_force as newton_second_law


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(1 * units.kilogram)
    a = QuantityVector([Quantity(3 * units.meter / units.second**2)])
    Args = namedtuple("Args", ["m", "a"])
    return Args(m=m, a=a)


def test_basic_force(test_args):
    result = newton_second_law.calculate_force(test_args.m, test_args.a)
    assert len(result.components) == 1
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force = convert_to(result.components[0], units.newton).evalf(2)
    assert result_force == approx(3.0, 0.01)
