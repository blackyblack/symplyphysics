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
from symplyphysics.core.vectors.vectors import QuantityVector
from symplyphysics.laws.dynamics.vector import instantaneous_power_is_force_dot_velocity as power_law

# Description
## A force is acting on an object, and at some time, the force vector is (1, 1, -1) N, and the object
## moves with the velocity vector equal to (2, 0, -1) m/s at that time. The power of the force exerted
## on the object is 3 W.


@fixture(name="test_args")
def test_args_fixture():
    f = QuantityVector([
        Quantity(1.0 * units.newton),
        Quantity(1.0 * units.newton),
        Quantity(-1.0 * units.newton),
    ])
    v = QuantityVector([
        Quantity(2.0 * units.meter / units.second),
        Quantity(0.0 * units.meter / units.second),
        Quantity(-1.0 * units.meter / units.second),
    ])
    Args = namedtuple("Args", "f v")
    return Args(f=f, v=v)


def test_basic_law(test_args):
    result = power_law.calculate_power(test_args.f, test_args.v)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.power)
    result_value = convert_to(result, units.watt).evalf(3)
    assert_approx(result_value, 3)


def test_bad_force(test_args):
    f_bad_vector = QuantityVector([
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
    ])
    with raises(errors.UnitsError):
        power_law.calculate_power(f_bad_vector, test_args.v)

    f_scalar = Quantity(1.0 * units.newton)
    with raises(AttributeError):
        power_law.calculate_power(f_scalar, test_args.v)

    with raises(TypeError):
        power_law.calculate_power(100, test_args.v)
    with raises(TypeError):
        power_law.calculate_power([100], test_args.v)


def test_bad_velocity(test_args):
    v_bad_vector = QuantityVector([
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
    ])
    with raises(errors.UnitsError):
        power_law.calculate_power(test_args.f, v_bad_vector)

    v_scalar = Quantity(1.0 * units.meter / units.second)
    with raises(AttributeError):
        power_law.calculate_power(test_args.f, v_scalar)

    with raises(TypeError):
        power_law.calculate_power(test_args.f, 100)
    with raises(TypeError):
        power_law.calculate_power(test_args.f, [100])
