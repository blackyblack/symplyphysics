from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.kinematics.vector import (
    displacement_is_angular_displacement_cross_radius as linear_displacement_law,)

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

# Description
## A body is rotating about a fixes axis. It makes a rotation of 1e-5 rad in the positive direction around the
## z-axis in the xy-plane, small enough so that the body's radius vector can be considered constant and equal to
## (0, 0.1, 0) m. During this rotation the body's linear displacement amounts to (-1e-6, 0, 0) m. Note that the angular
## displacement is a pseudovector, the magnitude of which is the angle of rotation and which is aligned along
## the axis of rotation. Its direction can be found via the right-hand rule.

Args = namedtuple("Args", "theta r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    theta = QuantityCoordinateVector([1e-5, 2e-5, -3e-5], CARTESIAN)

    r = QuantityCoordinateVector([
        Quantity(0.428571428571428 * units.meter),
        Quantity(0.857142857142857 * units.meter),
        Quantity(0.714285714285714 * units.meter),
    ], CARTESIAN)
    return Args(theta=theta, r=r)


def test_basic_law(test_args: Args) -> None:
    result = linear_displacement_law.calculate_linear_displacement(test_args.theta, test_args.r)

    expected = QuantityCoordinateVector([
        4e-5 * units.meter,
        -2e-5 * units.meter,
        0,
    ], CARTESIAN)

    assert_equal_vectors(result, expected, absolute_tolerance=1e-10)


def test_bad_angular_displacement(test_args: Args) -> None:
    theta_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        linear_displacement_law.calculate_linear_displacement(theta_bad_vector, test_args.r)

    theta_non_orthogonal = QuantityCoordinateVector([0, 1e-5, 1e-5], CARTESIAN)
    with raises(ValueError):
        linear_displacement_law.calculate_linear_displacement(theta_non_orthogonal, test_args.r)

    theta_scalar = Quantity(1.0 * units.radian)
    with raises(ValueError):
        linear_displacement_law.calculate_linear_displacement(theta_scalar, test_args.r)
    with raises(ValueError):
        linear_displacement_law.calculate_linear_displacement(100, test_args.r)
    with raises(ValueError):
        linear_displacement_law.calculate_linear_displacement([100], test_args.r)


def test_bad_rotation_radius(test_args: Args) -> None:
    r_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        linear_displacement_law.calculate_linear_displacement(test_args.theta, r_bad_vector)

    r_scalar = Quantity(1.0 * units.meter)
    with raises(ValueError):
        linear_displacement_law.calculate_linear_displacement(test_args.theta, r_scalar)

    with raises(TypeError):
        linear_displacement_law.calculate_linear_displacement(test_args.theta, 100)
    with raises(TypeError):
        linear_displacement_law.calculate_linear_displacement(test_args.theta, [100])
