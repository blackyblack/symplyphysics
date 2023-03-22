from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.definitions import wavelength_is_gap_between_peaks as wavelength_definition

# Description.
## Speed of light in air is 299910 km/s. Frequency of radio Europa+ is 101.6 MHz and therefore has oscillation period of 9,84252e-9 seconds. 
## According to online calculator https://vpayaem.ru/inf_wave1.html wavelength should be 2.95m.

@fixture
def test_args():
    v1 = units.Quantity('v1')
    SI.set_quantity_dimension(v1, units.velocity)
    SI.set_quantity_scale_factor(v1, 299910000 * units.meter / units.second)
    
    period1 = units.Quantity('period1')
    SI.set_quantity_dimension(period1, units.time)
    SI.set_quantity_scale_factor(period1, 9.84252e-9 * units.second)

    Args = namedtuple('Args', ['v1', 'period1'])
    return Args(v1=v1, period1=period1)

def test_basic_wavelength(test_args):
    result = wavelength_definition.calculate_wavelength(test_args.v1, test_args.period1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_wavelength = convert_to(result, wavelength_definition.definition_dimension_SI).subs({units.meter: 1}).evalf(6)    
    assert result_wavelength == approx(2.95, 0.001)

def test_bad_velocity(test_args):
    vb = units.Quantity('vb')
    SI.set_quantity_dimension(vb, units.length)
    SI.set_quantity_scale_factor(vb, 1 * units.meter)
    with raises(errors.UnitsError):
        wavelength_definition.calculate_wavelength(vb, test_args.period1)

    with raises(TypeError):
        wavelength_definition.calculate_wavelength(100, test_args.period1)

def test_bad_period(test_args):
    pb = units.Quantity('pb')
    SI.set_quantity_dimension(pb, units.length)
    SI.set_quantity_scale_factor(pb, 1 * units.meter)

    with raises(errors.UnitsError):
        wavelength_definition.calculate_wavelength(test_args.v1, pb)

    with raises(TypeError):
        wavelength_definition.calculate_wavelength(test_args.v1, 100)
