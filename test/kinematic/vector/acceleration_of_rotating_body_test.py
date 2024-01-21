from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
    QuantityVector,
)
from symplyphysics.laws.kinematic.vector import acceleration_of_rotating_body as acceleration_law

# Description
## A body is moving along a curve. At a certain point in time, its radial acceleration
## with respect to the momentary rotational axis is (-1, 0, 0.4) m/s^2, and its tangential
## acceleration is (0, 0.1, -3) m/s^2. Its acceleration should amount to (-1.0, 0.1, -2.6) m/s^2.


@fixture(name="test_args")
def test_args_fixture():
    a_r = QuantityVector([
        Quantity(-1.0 * units.meter / units.second**2),
        Quantity(0.0 * units.meter / units.second**2),
        Quantity(0.4 * units.meter / units.second**2),
    ])
    a_t = QuantityVector([
        Quantity(0.0 * units.meter / units.second**2),
        Quantity(0.1 * units.meter / units.second**2),
        Quantity(-3.0 * units.meter / units.second**2),
    ])
    Args = namedtuple("Args", "a_r a_t")
    return Args(a_r=a_r, a_t=a_t)


def test_basic_law(test_args):
    result = acceleration_law.calculate_acceleration(test_args.a_r, test_args.a_t)
    assert len(result.components) == 3
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.acceleration)
    correct_values = (-1.0, 0.1, -2.6)
    for result_component, correct_value in zip(result.components, correct_values):
        result_value = convert_to(result_component, units.meter / units.second**2).evalf(3)
        assert result_value == approx(correct_value, 1e-3)


def test_bad_acceleration(test_args):
    a_bad_vector = QuantityVector([
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
    ])
    with raises(errors.UnitsError):
        acceleration_law.calculate_acceleration(a_bad_vector, test_args.a_t)
    with raises(errors.UnitsError):
        acceleration_law.calculate_acceleration(test_args.a_r, a_bad_vector)

    a_bad_scalar = Quantity(1.0 * units.meter / units.second**2)
    with raises(AttributeError):
        acceleration_law.calculate_acceleration(a_bad_scalar, test_args.a_t)
    with raises(AttributeError):
        acceleration_law.calculate_acceleration(test_args.a_r, a_bad_scalar)

    with raises(TypeError):
        acceleration_law.calculate_acceleration(100, test_args.a_t)
    with raises(TypeError):
        acceleration_law.calculate_acceleration(test_args.a_r, 100)
