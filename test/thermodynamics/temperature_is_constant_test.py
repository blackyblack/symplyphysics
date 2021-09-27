from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.thermodynamics import temperature_is_constant as boyles_law

@fixture
def test_args():
    P0 = units.Quantity('P0')
    SI.set_quantity_dimension(P0, units.pressure)
    SI.set_quantity_scale_factor(P0, 1 * units.pascal)
    P1 = units.Quantity('P1')
    SI.set_quantity_dimension(P1, units.pressure)
    SI.set_quantity_scale_factor(P1, 2 * units.pascal)
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.volume)
    SI.set_quantity_scale_factor(V0, 1 * units.liter)

    Args = namedtuple('Args', ['P0', 'P1', 'V0'])
    return Args(P0 = P0, P1 = P1, V0 = V0)

def test_basic_volume(test_args):
    result = boyles_law.calculate_volume(test_args.P0, test_args.V0, test_args.P1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.volume)

    result_volume = convert_to(result, units.liter).subs(units.liter, 1).evalf(2)
    assert result_volume == approx(0.5, 0.01)

def test_bad_pressure(test_args):
    P0b = units.Quantity('P0b')
    SI.set_quantity_dimension(P0b, units.length)
    SI.set_quantity_scale_factor(P0b, 1 * units.meter)
    with raises(errors.UnitsError):
        boyles_law.calculate_volume(P0b, test_args.V0, test_args.P1)

    # Make P1 invalid
    P1b = units.Quantity('P1b')
    SI.set_quantity_dimension(P1b, units.length)
    SI.set_quantity_scale_factor(P1b, 1 * units.meter)
    with raises(errors.UnitsError):
        boyles_law.calculate_volume(test_args.P0, test_args.V0, P1b)

    with raises(TypeError):
        boyles_law.calculate_volume(100, test_args.V0, test_args.P1)
    
    with raises(TypeError):
        boyles_law.calculate_volume(test_args.P0, test_args.V0, 100)

def test_bad_volume(test_args):
    V0b = units.Quantity('V0b')
    SI.set_quantity_dimension(V0b, units.length)
    SI.set_quantity_scale_factor(V0b, 1 * units.meter)

    with raises(errors.UnitsError):
        boyles_law.calculate_volume(test_args.P0, V0b, test_args.P1)

    with raises(TypeError):
        boyles_law.calculate_volume(test_args.P0, 100, test_args.P1)
