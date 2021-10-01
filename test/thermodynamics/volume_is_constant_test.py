from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.thermodynamics import volume_is_constant as isochoric_law

@fixture
def test_args():
    t0 = units.Quantity('t0')
    SI.set_quantity_dimension(t0, units.temperature)
    SI.set_quantity_scale_factor(t0, 1 * units.kelvin)
    t1 = units.Quantity('t1')
    SI.set_quantity_dimension(t1, units.temperature)
    SI.set_quantity_scale_factor(t1, 2 * units.kelvin)
    P0 = units.Quantity('P0')
    SI.set_quantity_dimension(P0, units.pressure)
    SI.set_quantity_scale_factor(P0, 1 * units.pascal)

    Args = namedtuple('Args', ['t0', 't1', 'P0'])
    return Args(t0 = t0, t1 = t1, P0 = P0)

def test_basic_pressure(test_args):
    result = isochoric_law.calculate_pressure(test_args.t0, test_args.P0, test_args.t1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)

    result_pressure = convert_to(result, units.pascal).subs(units.pascal, 1).evalf(2)
    assert result_pressure == approx(2.0, 0.01)

def test_bad_temperature(test_args):
    # Make t0 invalid
    t0b = units.Quantity('t0b')
    SI.set_quantity_dimension(t0b, units.length)
    SI.set_quantity_scale_factor(t0b, 1 * units.meter)
    with raises(errors.UnitsError):
        isochoric_law.calculate_pressure(t0b, test_args.P0, test_args.t1)

    # Make t1 invalid
    t1b = units.Quantity('t1b')
    SI.set_quantity_dimension(t1b, units.length)
    SI.set_quantity_scale_factor(t1b, 1 * units.meter)
    with raises(errors.UnitsError):
        isochoric_law.calculate_pressure(test_args.t0, test_args.P0, t1b)

    with raises(TypeError):
        isochoric_law.calculate_pressure(100, test_args.P0, test_args.t1)
    
    with raises(TypeError):
        isochoric_law.calculate_pressure(test_args.t0, test_args.P0, 100)

def test_bad_pressure(test_args):
    P0b = units.Quantity('P0b')
    SI.set_quantity_dimension(P0b, units.length)
    SI.set_quantity_scale_factor(P0b, 1 * units.meter)

    with raises(errors.UnitsError):
        isochoric_law.calculate_pressure(test_args.t0, P0b, test_args.t1)

    with raises(TypeError):
        isochoric_law.calculate_pressure(test_args.t0, 100, test_args.t1)
