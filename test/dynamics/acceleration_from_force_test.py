from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import acceleration_from_force as newton_second_law


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(1 * units.kilogram)
    a = Quantity(3 * units.meter / units.second**2)
    Args = namedtuple("Args", ["m", "a"])
    return Args(m=m, a=a)


def test_basic_force(test_args):
    result = newton_second_law.calculate_force(test_args.m, test_args.a)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force = convert_to(result, units.newton).evalf(2)
    assert_approx(result_force, 3)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_force(mb, test_args.a)
    with raises(TypeError):
        newton_second_law.calculate_force(100, test_args.a)


def test_bad_acceleration(test_args):
    ab = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_force(test_args.m, ab)
    with raises(TypeError):
        newton_second_law.calculate_force(test_args.m, 100)
