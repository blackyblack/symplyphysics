from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration

@fixture
def test_args():
    v0 = units.Quantity('v0')
    SI.set_quantity_dimension(v0, units.velocity)
    SI.set_quantity_scale_factor(v0, 1 * units.meter / units.second)
    v1 = units.Quantity('v1')
    SI.set_quantity_dimension(v1, units.velocity)
    SI.set_quantity_scale_factor(v1, 20 * units.meter / units.second)
    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.time)
    SI.set_quantity_scale_factor(t, 5 * units.second)

    Args = namedtuple('Args', ['v0', 'v1', 't'])
    return Args(v0 = v0, v1 = v1, t = t)

def test_basic_acceleration(test_args):
    result = acceleration.calculate_linear_acceleration(test_args.v0, test_args.v1, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.acceleration)

    result_acceleration = convert_to(result, acceleration.definition_dimension_SI).subs({
        units.meter: 1, units.second: 1}).evalf(2)
    assert result_acceleration == approx(3.8, 0.01)

def test_bad_velocity(test_args):
    v0b = units.Quantity('v0b')
    SI.set_quantity_dimension(v0b, units.length)
    SI.set_quantity_scale_factor(v0b, 1 * units.meter)
    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(v0b, test_args.v1, test_args.t)

    # Make v1 invalid
    v1b = units.Quantity('v1b')
    SI.set_quantity_dimension(v1b, units.length)
    SI.set_quantity_scale_factor(v1b, 1 * units.meter)
    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(test_args.v0, v1b, test_args.t)

    with raises(TypeError):
        acceleration.calculate_linear_acceleration(100, test_args.v1, test_args.t)

    with raises(TypeError):
        acceleration.calculate_linear_acceleration(test_args.v0, 100, test_args.t)

def test_bad_time(test_args):
    tb = units.Quantity('tb')
    SI.set_quantity_dimension(tb, units.length)
    SI.set_quantity_scale_factor(tb, 1 * units.meter)

    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(test_args.v0, test_args.v1, tb)
    
    with raises(TypeError):
        acceleration.calculate_linear_acceleration(test_args.v0, test_args.v1, 100)
