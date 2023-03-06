from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.definitions import wavelength_is_velocity_by_frequency as wavelength_definition

# Description.
## Spreading speed of light is the air is 299704 m/s. Refraction factor of air is 1.003. Frequency of radio Europa+ is 101.6 MHz. 
## According to online calculator https://vpayaem.ru/inf_wave1.html wavelength should be 2.94mm.

@fixture
def test_args():
    v = units.Quantity('v')
    SI.set_quantity_dimension(v, units.velocity)
    SI.set_quantity_scale_factor(v, 299704 * units.meter / units.second)
    
    refraction_factor = 1.003

    frequency = units.Quantity('frequency')
    SI.set_quantity_dimension(frequency, units.frequency)
    SI.set_quantity_scale_factor(frequency, 101.6 * 1000000 * units.hertz)

    Args = namedtuple('Args', ['v', 'refraction_factor', 'frequency'])
    return Args(v=v, refraction_factor=refraction_factor, frequency=frequency)

def test_basic_wavelength(test_args):
    result = wavelength_definition.calculate_wavelength(test_args.v, test_args.refraction_factor, test_args.frequency)
    assert SI.get_dimension_system().equivalent_dims(
        result.dimension, units.length)

    result_wavelength = convert_to(result, wavelength_definition.definition_dimension_SI).subs({
        units.meter: 1}).evalf(3)
    assert result_wavelength == approx(0.00294, 0.001)


def test_bad_velocity(test_args):
    vb = units.Quantity('vb')
    SI.set_quantity_dimension(vb, units.length)
    SI.set_quantity_scale_factor(vb, 1 * units.meter)
    with raises(errors.UnitsError):
        wavelength_definition.calculate_wavelength(vb, test_args.refraction_factor, test_args.frequency)

    with raises(TypeError):
        wavelength_definition.calculate_wavelength(100, test_args.refraction_factor, test_args.frequency)


def test_bad_frequency(test_args):
    fb = units.Quantity('fb')
    SI.set_quantity_dimension(fb, units.length)
    SI.set_quantity_scale_factor(fb, 1 * units.meter)

    with raises(errors.UnitsError):
        wavelength_definition.calculate_wavelength(test_args.v, test_args.refraction_factor, fb)

    with raises(TypeError):
        wavelength_definition.calculate_wavelength(test_args.v, test_args.refraction_factor, 100)

