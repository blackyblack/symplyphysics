from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
    dimensionless,
)
from symplyphysics.laws.kinematic import (
    constant_angular_acceleration_rotation_is_parabolic as angular_displacement_law,)

# Description
## A body is rotating about a fixed axis with constant angular acceleration of -2 rad/s**2.
## Its initial angular velocity was 1 rad/s. In 5 s its angular displacement amounts to
## -20 rad.


@fixture(name="test_args")
def test_args_fixture():
    w0 = Quantity(1.0 * units.radian / units.second)
    alpha = Quantity(-2.0 * units.radian / units.second**2)
    t = Quantity(5.0 * units.second)
    Args = namedtuple("Args", "w0 alpha t")
    return Args(w0=w0, alpha=alpha, t=t)


def test_basic_law(test_args):
    result = angular_displacement_law.calculate_angular_displacement(test_args.w0, test_args.alpha,
        test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, dimensionless)
    result_value = convert_to(result, units.radian).evalf(3)
    assert result_value == approx(-20.0, 1e-3)


def test_bad_angular_velocity(test_args):
    w0b = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_displacement_law.calculate_angular_displacement(w0b, test_args.alpha, test_args.t)
    with raises(TypeError):
        angular_displacement_law.calculate_angular_displacement(100, test_args.alpha, test_args.t)


def test_bad_angular_acceleration(test_args):
    alpha_b = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_displacement_law.calculate_angular_displacement(test_args.w0, alpha_b, test_args.t)
    with raises(TypeError):
        angular_displacement_law.calculate_angular_displacement(test_args.w0, 100, test_args.t)


def test_bad_time(test_args):
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_displacement_law.calculate_angular_displacement(test_args.w0, test_args.alpha, tb)
    with raises(TypeError):
        angular_displacement_law.calculate_angular_displacement(test_args.w0, test_args.alpha, 100)
