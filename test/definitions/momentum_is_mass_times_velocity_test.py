from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.definitions import momentum_is_mass_times_velocity as momentum_def

@fixture
def test_args():
    m = units.Quantity('m')
    SI.set_quantity_dimension(m, units.mass)
    SI.set_quantity_scale_factor(m, 1 * units.kilogram)
    v = units.Quantity('v')
    SI.set_quantity_dimension(v, units.velocity)
    SI.set_quantity_scale_factor(v, 5 * units.meter / units.second)

    Args = namedtuple('Args', ['m', 'v'])
    return Args(m = m, v = v)

def test_basic_momentum(test_args):
    result = momentum_def.calculate_momentum(
    test_args.m, test_args.v)
    assert SI.get_dimension_system().equivalent_dims(
        result.dimension, units.momentum)
    result_momentum = convert_to(result, momentum_def.definition_dimension_SI).subs({
        units.kilogram * units.meter / units.second: 1}).evalf(2)    
    assert result_momentum == approx(5.0, 0.001)

def test_bad_mass(test_args):
    mb = units.Quantity('mb')
    SI.set_quantity_dimension(mb, units.length)
    SI.set_quantity_scale_factor(mb, 1 * units.meter)

    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(mb, test_args.v)

    with raises(TypeError):
        momentum_def.calculate_momentum(100, test_args.v)

def test_bad_velocity(test_args):
    vb = units.Quantity('vb')
    SI.set_quantity_dimension(vb, units.length)
    SI.set_quantity_scale_factor(vb, 3 * units.meter)

    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(test_args.m, vb)

    with raises(TypeError):
        momentum_def.calculate_momentum(test_args.m, 100)
