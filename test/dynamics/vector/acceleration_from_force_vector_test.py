from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    units,
    errors,
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
    f = QuantityVector([Quantity(3 * units.newton)])
    Args = namedtuple("Args", ["m", "a", "f"])
    return Args(m=m, a=a, f=f)


def test_basic_force(test_args):
    result = newton_second_law.calculate_force(test_args.m, test_args.a)
    assert len(result.components) == 1
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force = convert_to(result.components[0], units.newton).evalf(2)
    assert_approx(result_force, 3)


def test_basic_acceleration(test_args):
    result = newton_second_law.calculate_acceleration(test_args.m, test_args.f)
    assert len(result.components) == 1
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.acceleration)
    result_acceleration = convert_to(result.components[0], units.meter / units.second**2).evalf(2)
    assert_approx(result_acceleration, 3)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_force(mb, test_args.a)
    with raises(TypeError):
        newton_second_law.calculate_force(100, test_args.a)
    with raises(errors.UnitsError):
        newton_second_law.calculate_acceleration(mb, test_args.f)
    with raises(TypeError):
        newton_second_law.calculate_acceleration(100, test_args.f)


def test_bad_acceleration(test_args):
    ab = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_force(test_args.m, ab)
    with raises(TypeError):
        newton_second_law.calculate_force(test_args.m, 100)


def test_bad_force(test_args):
    fb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_acceleration(test_args.m, fb)
    with raises(TypeError):
        newton_second_law.calculate_acceleration(test_args.m, 100)
