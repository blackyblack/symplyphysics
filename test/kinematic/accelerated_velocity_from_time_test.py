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
    initial_velocity1 = units.Quantity('initial_velocity1')
    SI.set_quantity_dimension(initial_velocity1, units.velocity)
    SI.set_quantity_scale_factor(initial_velocity1, 2 * units.meter / units.second)

    acceleration = units.Quantity('acceleration')
    SI.set_quantity_dimension(acceleration, units.acceleration)
    SI.set_quantity_scale_factor(acceleration, -9.8 * units.meter / units.second**2)

    time = units.Quantity('time')
    SI.set_quantity_dimension(time, units.time)
    SI.set_quantity_scale_factor(time, 5 * units.second)

    Args = namedtuple('Args', ['initial_velocity1', 'acceleration', 'time'])
    return Args(initial_velocity1 = initial_velocity1, acceleration = acceleration, time = time)

def test_basic_velocity(test_args):
    result = accelerated_velocity_law.calculate_velocity(test_args.initial_velocity1, test_args.acceleration, test_args.time)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity)
    result_vector = convert_to(result, units.meter/units.second).subs(units.meter/units.second, 1).evalf(2)
    assert result_vector == approx(-47, 0.01)

'''
def test_bad_angle(test_args):    
    ab = units.Quantity('Rb')
    SI.set_quantity_dimension(ab, units.length)
    SI.set_quantity_scale_factor(ab, 1 * units.meter)

    with raises(errors.UnitsError):
        projection_law.calculate_projection(test_args.force_vector_amplitude, ab)

    with raises(TypeError):
        projection_law.calculate_projection(test_args.force_vector_amplitude, 100)

        '''