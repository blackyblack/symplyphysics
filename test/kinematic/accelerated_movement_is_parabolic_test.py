from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)

from symplyphysics.laws.kinematic import accelerated_movement_is_parabolic as movement_law
# Description
## The man starts running with initial speed 10m/s and is getting tired while running (decreases his velocity) with -0.1 m/s/s. How far this man can run away in 10 seconds?

@fixture
def test_args():
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.velocity)
    SI.set_quantity_scale_factor(V0, 10 * units.meter / units.second)

    a = units.Quantity('a')
    SI.set_quantity_dimension(a, units.acceleration)
    SI.set_quantity_scale_factor(a, -0.1 * units.meter / units.second**2)

    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.time)
    SI.set_quantity_scale_factor(t, 10 * units.second)

    Args = namedtuple('Args', ['V0', 'a', 't'])
    return Args(V0 = V0, a = a, t = t)

def test_basic_distance(test_args):
    result = movement_law.calculate_distance(test_args.V0, test_args.a, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_vector = convert_to(result, units.meter).subs(units.meter, 1).evalf(2)
    assert result_vector == approx(95, 0.01)

'''
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
        '''