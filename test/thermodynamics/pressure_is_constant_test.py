from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.thermodynamics import pressure_is_constant as gay_lussacs_law

@fixture
def test_args():
    t0 = units.Quantity('t0')
    SI.set_quantity_dimension(t0, units.temperature)
    SI.set_quantity_scale_factor(t0, 1 * units.kelvin)
    t1 = units.Quantity('t1')
    SI.set_quantity_dimension(t1, units.temperature)
    SI.set_quantity_scale_factor(t1, 2 * units.kelvin)
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.volume)
    SI.set_quantity_scale_factor(V0, 1 * units.liter)

    Args = namedtuple('Args', ['t0', 't1', 'V0'])
    return Args(t0 = t0, t1 = t1, V0 = V0)

def test_basic_volume(test_args):
    result = gay_lussacs_law.calculate_volume(test_args.t0, test_args.V0, test_args.t1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.volume)

    result_volume = convert_to(result, units.liter).subs(units.liter, 1).evalf(2)
    assert result_volume == approx(2.0, 0.01)

def test_bad_temperature(test_args):
    # Make t0 invalid
    t0b = units.Quantity('t0b')
    SI.set_quantity_dimension(t0b, units.length)
    SI.set_quantity_scale_factor(t0b, 1 * units.meter)
    with raises(errors.UnitsError):
        gay_lussacs_law.calculate_volume(t0b, test_args.V0, test_args.t1)

    # Make t1 invalid
    t1b = units.Quantity('t1b')
    SI.set_quantity_dimension(t1b, units.length)
    SI.set_quantity_scale_factor(t1b, 1 * units.meter)
    with raises(errors.UnitsError):
        gay_lussacs_law.calculate_volume(test_args.t0, test_args.V0, t1b)

    with raises(TypeError):
        gay_lussacs_law.calculate_volume(100, test_args.V0, test_args.t1)
    
    with raises(TypeError):
        gay_lussacs_law.calculate_volume(test_args.t0, test_args.V0, 100)

def test_bad_volume(test_args):
    V0b = units.Quantity('V0b')
    SI.set_quantity_dimension(V0b, units.length)
    SI.set_quantity_scale_factor(V0b, 1 * units.meter)

    with raises(errors.UnitsError):
        gay_lussacs_law.calculate_volume(test_args.t0, V0b, test_args.t1)

    with raises(TypeError):
        gay_lussacs_law.calculate_volume(test_args.t0, 100, test_args.t1)
