from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    SI,
    QuantityVector,
)
from symplyphysics.laws.kinematic.vector import (
    rolling_velocity_is_translational_plus_rotational_velocities as rolling_velocity_law,
)

# Description
## A wheel is rolling smoothly, and for a point on the wheel the translational velocity
## component amounts to (1, 2, -1) m/s, and the rotational one to (0, -1, 2) m/s.
## Then the linear velocity of that point is (1, 1, 1) m/s at that point in time.


@fixture(name="test_args")
def test_args_fixture():
    v_t = QuantityVector([
        Quantity(1.0 * units.meter / units.second),
        Quantity(2.0 * units.meter / units.second),
        Quantity(-1.0 * units.meter / units.second),
    ])
    v_r = QuantityVector([
        Quantity(0.0 * units.meter / units.second),
        Quantity(-1.0 * units.meter / units.second),
        Quantity(2.0 * units.meter / units.second),
    ])
    Args = namedtuple("Args", "v_t v_r")
    return Args(v_t=v_t, v_r=v_r)


def test_law(test_args):
    result = rolling_velocity_law.calculate_rolling_velocity(test_args.v_t, test_args.v_r)
    assert len(result.components) == 3
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity)
    for result_component, correct_value in zip(result.components, [1.0, 1.0, 1.0]):
        assert_equal(result_component, correct_value * units.meter / units.second)


def test_bad_velocities(test_args):
    v_bad_vector = QuantityVector([
        Quantity(1.0 * units.meter),
        Quantity(2.0 * units.meter),
        Quantity(-1.0 * units.meter),
    ])
    with raises(errors.UnitsError):
        rolling_velocity_law.calculate_rolling_velocity(v_bad_vector, test_args.v_r)
    with raises(errors.UnitsError):
        rolling_velocity_law.calculate_rolling_velocity(test_args.v_t, v_bad_vector)

    v_scalar = Quantity(1.0 * units.meter / units.second)
    with raises(AttributeError):
        rolling_velocity_law.calculate_rolling_velocity(v_scalar, test_args.v_r)
    with raises(AttributeError):
        rolling_velocity_law.calculate_rolling_velocity(test_args.v_t, v_scalar)

    with raises(TypeError):
        rolling_velocity_law.calculate_rolling_velocity(100, test_args.v_r)
    with raises(TypeError):
        rolling_velocity_law.calculate_rolling_velocity([100], test_args.v_r)
    with raises(TypeError):
        rolling_velocity_law.calculate_rolling_velocity(test_args.v_t, 100)
    with raises(TypeError):
        rolling_velocity_law.calculate_rolling_velocity(test_args.v_t, [100])
