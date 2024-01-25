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
    linear_displacement_is_angular_displacement_cross_radius as linear_displacement_law,
)

# Description
## A body is rotating about a fixes axis. It makes a rotation of 1e-5 rad in the positive direction around the
## z-axis in the xy-plane, small enough so that the body's radius vector can be considered constant and equal to 
## (0, 0.1, 0) m. During this rotation the body's linear displacement amounts to (-1e-6, 0, 0) m. Note that the angular
## displacement is a pseudovector, the magnitude of which is the angle of rotation and which is aligned along
## the axis of rotation. Its direction can be found via the right-hand rule.


@fixture(name="test_args")
def test_args_fixture():
    theta = QuantityVector([0.0, 0.0, 1e-5])
    r = QuantityVector([
        Quantity(0.0 * units.meter),
        Quantity(0.1 * units.meter),
        Quantity(0.0 * units.meter),
    ])
    Args = namedtuple("Args", "theta r")
    return Args(theta=theta, r=r)


def test_basic_law(test_args):
    result = linear_displacement_law.calculate_linear_displacement(test_args.theta, test_args.r)
    assert len(result.components) == 3
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    for result_component, correct_value in zip(result.components, [-1e-6, 0.0, 0.0]):
        result_value = convert_to(result_component, units.meter).evalf(3)
        assert result_value == approx(correct_value, 1e-4)


def test_bad_angular_displacement(test_args):
    theta_bad_vector = QuantityVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ])
    with raises(errors.UnitsError):
        linear_displacement_law.calculate_linear_displacement(theta_bad_vector, test_args.r)

    theta_non_orthogonal = QuantityVector([0.0, 1e-5, 1e-5])
    with raises(ValueError):
        linear_displacement_law.calculate_linear_displacement(theta_non_orthogonal, test_args.r)

    theta_scalar = Quantity(1.0 * units.radian)
    with raises(AttributeError):
        linear_displacement_law.calculate_linear_displacement(theta_scalar, test_args.r)
    with raises(AttributeError):
        linear_displacement_law.calculate_linear_displacement(100, test_args.r)
    with raises(AttributeError):
        linear_displacement_law.calculate_linear_displacement([100], test_args.r)


def test_bad_rotation_radius(test_args):
    r_bad_vector = QuantityVector([
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
    ])
    with raises(errors.UnitsError):
        linear_displacement_law.calculate_linear_displacement(test_args.theta, r_bad_vector)

    r_scalar = Quantity(1.0 * units.meter)
    with raises(AttributeError):
        linear_displacement_law.calculate_linear_displacement(test_args.theta, r_scalar)

    with raises(TypeError):
        linear_displacement_law.calculate_linear_displacement(test_args.theta, 100)
    with raises(TypeError):
        linear_displacement_law.calculate_linear_displacement(test_args.theta, [100])
