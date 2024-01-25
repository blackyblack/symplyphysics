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
from symplyphysics.laws.kinematic.vector import (
    linear_velocity_is_angular_velocity_cross_radius as linear_velocity_law,
)

# Description
## A rigid body is moving about an axis with angular velocity 4.0 rad/s in the positive direction
## of the z-axis (therefore, the angular velocity pseudovector is (0.0, 0.0, 4.0) rad/s). At a
## certain point in time, the radius vector of the body is (0, -0.5, 0) m. Its linear velocity at
## that time is (2.0, 0, 0) m/s.


@fixture(name="test_args")
def test_args_fixture():
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
    Args = namedtuple("Args", "w r")
    return Args(w=w, r=r)


def test_basic_law(test_args):
    result = linear_velocity_law.calculate_linear_velocity(test_args.w, test_args.r)
    assert len(result.components) == 3
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity)
    for result_component, correct_value in zip(result.components, [2.0, 0.0, 0.0]):
        result_value = convert_to(result_component, units.meter / units.second).evalf(3)
        assert result_value == approx(correct_value, 1e-3)


def test_bad_angular_velocity(test_args):
    w_bad_vector = QuantityVector([
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
    ])
    with raises(errors.UnitsError):
        linear_velocity_law.calculate_linear_velocity(w_bad_vector, test_args.r)

    w_scalar = Quantity(1.0 * units.radian / units.second)
    with raises(AttributeError):
        linear_velocity_law.calculate_linear_velocity(w_scalar, test_args.r)

    with raises(TypeError):
        linear_velocity_law.calculate_linear_velocity(100, test_args.r)
    with raises(TypeError):
        linear_velocity_law.calculate_linear_velocity([100], test_args.r)


def test_bad_rotation_radius(test_args):
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
