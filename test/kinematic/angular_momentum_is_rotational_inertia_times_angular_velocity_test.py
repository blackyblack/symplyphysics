from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic import (
    angular_momentum_is_rotational_inertia_times_angular_velocity as angular_momentum_law,)

# Description
## A rigid body with rotational inertia of 1.5 kg*m**2 is rotating about a fixed axis
## with angular velocity of 4.0 rad/s. Its angular momentum amounts to 6.0 kg*m**2/s.


@fixture(name="test_args")
def test_args_fixture():
    I = Quantity(1.5 * units.kilogram * units.meter**2)
    w = Quantity(4.0 * units.radian / units.second)
    Args = namedtuple("Args", "I w")
    return Args(I=I, w=w)


def test_law(test_args):
    result = angular_momentum_law.calculate_angular_momentum(test_args.I, test_args.w)
    assert_equal(result, 6.0 * units.kilogram * units.meter**2 / units.second)


def test_bad_rotational_inertia(test_args):
    Ib = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        angular_momentum_law.calculate_angular_momentum(Ib, test_args.w)
    with raises(TypeError):
        angular_momentum_law.calculate_angular_momentum(100, test_args.w)


def test_bad_angular_velocity(test_args):
    wb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        angular_momentum_law.calculate_angular_momentum(test_args.I, wb)
    with raises(TypeError):
        angular_momentum_law.calculate_angular_momentum(test_args.I, 100)
