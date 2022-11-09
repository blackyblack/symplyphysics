from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)

from symplyphysics.laws.kinematic import accelerated_velocity_from_time as accelerated_velocity_law
# Description
## We are having object falling with initial speed 2 m/s directed upwards and 9.8 m/s**2 gravitation acceleration.
## Calculate speed after 5 secs of flight.
## Let's choose Y axis directed upwards and space is 1-dimensional. So initial speed is +2m/s (as the angle between velocity and Y is 0) and acceleration is -9.8 m/s**2 (as the angle between acceleration and Y is 180)

@fixture
def test_args():
    V1 = units.Quantity('V1')
    SI.set_quantity_dimension(V1, units.velocity)
    SI.set_quantity_scale_factor(V1, 2 * units.meter / units.second)

    A1 = units.Quantity('A1')
    SI.set_quantity_dimension(A1, units.acceleration)
    SI.set_quantity_scale_factor(A1, -9.8 * units.meter / units.second**2)

    T1 = units.Quantity('T1')
    SI.set_quantity_dimension(T1, units.time)
    SI.set_quantity_scale_factor(T1, 5 * units.second)

    Args = namedtuple('Args', ['V1', 'A1', 'T1'])
    return Args(V1 = V1, A1 = A1, T1 = T1)

def test_basic_velocity(test_args):
    result = accelerated_velocity_law.calculate_velocity(test_args.V1, test_args.A1, test_args.T1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity)
    result_vector = convert_to(result, units.meter/units.second).subs(units.meter/units.second, 1).evalf(2)
    assert result_vector == approx(-47, 0.01)

def test_bad_velocity(test_args):
    Vb = units.Quantity('Vb')
    SI.set_quantity_dimension(Vb, units.length)
    SI.set_quantity_scale_factor(Vb, 1 * units.meter)

    with raises(errors.UnitsError):
        accelerated_velocity_law.calculate_velocity(Vb, test_args.A1, test_args.T1)

    with raises(TypeError):
        accelerated_velocity_law.calculate_velocity(100, test_args.A1, test_args.T1)

def test_bad_acceleration(test_args):
    Ab = units.Quantity('Ab')
    SI.set_quantity_dimension(Ab, units.length)
    SI.set_quantity_scale_factor(Ab, 1 * units.meter)

    with raises(errors.UnitsError):
        accelerated_velocity_law.calculate_velocity(test_args.V1, Ab, test_args.T1)

    with raises(TypeError):
        accelerated_velocity_law.calculate_velocity(test_args.V1, 100, test_args.T1)   

def test_bad_time(test_args):
    Tb = units.Quantity('Tb')
    SI.set_quantity_dimension(Tb, units.length)
    SI.set_quantity_scale_factor(Tb, 1 * units.meter)

    with raises(errors.UnitsError):
        accelerated_velocity_law.calculate_velocity(test_args.V1, test_args.A1, Tb)

    with raises(TypeError):
        accelerated_velocity_law.calculate_velocity(test_args.V1, test_args.A1, 100)              