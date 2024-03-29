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
    linear_velocity_is_angular_velocity_cross_radius as linear_velocity_law,)

# Description
## A rigid body is moving about an axis with angular velocity 4.0 rad/s in the positive direction around
## the z-axis in the xy-plane (therefore, the angular velocity pseudovector is (0.0, 0.0, 4.0) rad/s). At a
## certain point in time, the radius vector of the body is (0, -0.5, 0) m. Its linear velocity at
## that time is (2.0, 0, 0) m/s.

Args = namedtuple("Args", "w r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w = QuantityVector([
        Quantity(0.0 * units.radian / units.second),
        Quantity(0.0 * units.radian / units.second),
        Quantity(4.0 * units.radian / units.second),
    ])
    r = QuantityVector([
        Quantity(0.0 * units.meter),
        Quantity(-0.5 * units.meter),
        Quantity(0.0 * units.meter),
    ])
    return Args(w=w, r=r)


def test_basic_law(test_args: Args) -> None:
    result = linear_velocity_law.calculate_linear_velocity(test_args.w, test_args.r)
    assert len(result.components) == 3
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity)
    for result_component, correct_value in zip(result.components, [2.0, 0.0, 0.0]):
        assert_equal(result_component, correct_value * units.meter / units.second)


def test_bad_angular_velocity(test_args: Args) -> None:
    w_bad_vector = QuantityVector([
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
    ])
    with raises(errors.UnitsError):
        linear_velocity_law.calculate_linear_velocity(w_bad_vector, test_args.r)

    w_non_orthogonal = QuantityVector([
        Quantity(0.0 * units.radian / units.second),
        Quantity(1.0 * units.radian / units.second),
        Quantity(1.0 * units.radian / units.second),
    ])
    with raises(ValueError):
        linear_velocity_law.calculate_linear_velocity(w_non_orthogonal, test_args.r)

    w_scalar = Quantity(1.0 * units.radian / units.second)
    with raises(AttributeError):
        linear_velocity_law.calculate_linear_velocity(w_scalar, test_args.r)

    with raises(TypeError):
        linear_velocity_law.calculate_linear_velocity(100, test_args.r)
    with raises(TypeError):
        linear_velocity_law.calculate_linear_velocity([100], test_args.r)


def test_bad_rotation_radius(test_args: Args) -> None:
    r_bad_vector = QuantityVector([
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
    ])
    with raises(errors.UnitsError):
        linear_velocity_law.calculate_linear_velocity(test_args.w, r_bad_vector)

    r_scalar = Quantity(1.0 * units.meter)
    with raises(AttributeError):
        linear_velocity_law.calculate_linear_velocity(test_args.w, r_scalar)

    with raises(TypeError):
        linear_velocity_law.calculate_linear_velocity(test_args.w, 100)
    with raises(TypeError):
        linear_velocity_law.calculate_linear_velocity(test_args.w, [100])
