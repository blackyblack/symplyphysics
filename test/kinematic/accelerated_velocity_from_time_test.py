from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import accelerated_velocity_from_time as accelerated_velocity_law
# Description
## We are having object falling with initial speed 2 m/s directed upwards and 9.8 m/s**2 gravitation acceleration.
## Calculate speed after 5 secs of flight.
## Let's choose Y axis directed upwards and space is 1-dimensional. So initial speed is +2m/s (as the angle between velocity and Y is 0) and acceleration is -9.8 m/s**2 (as the angle between acceleration and Y is 180)


@fixture
def test_args():
    V1 = Quantity(2 * units.meter / units.second)
    A1 = Quantity(-9.8 * units.meter / units.second**2)
    T1 = Quantity(5 * units.second)
    Args = namedtuple("Args", ["V1", "A1", "T1"])
    return Args(V1=V1, A1=A1, T1=T1)


def test_basic_velocity(test_args):
    result = accelerated_velocity_law.calculate_velocity(test_args.V1, test_args.A1, test_args.T1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity)
    result_velocity = convert_to(result,
        units.meter / units.second).subs(units.meter / units.second, 1).evalf(2)
    assert result_velocity == approx(-47, 0.01)


def test_bad_velocity(test_args):
    Vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        accelerated_velocity_law.calculate_velocity(Vb, test_args.A1, test_args.T1)
    with raises(TypeError):
        accelerated_velocity_law.calculate_velocity(100, test_args.A1, test_args.T1)


def test_bad_acceleration(test_args):
    Ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        accelerated_velocity_law.calculate_velocity(test_args.V1, Ab, test_args.T1)
    with raises(TypeError):
        accelerated_velocity_law.calculate_velocity(test_args.V1, 100, test_args.T1)


def test_bad_time(test_args):
    Tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        accelerated_velocity_law.calculate_velocity(test_args.V1, test_args.A1, Tb)
    with raises(TypeError):
        accelerated_velocity_law.calculate_velocity(test_args.V1, test_args.A1, 100)
