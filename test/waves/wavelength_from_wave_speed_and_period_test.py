from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors, units, Quantity, SI, convert_to,
)
from symplyphysics.laws.waves import wavelength_from_wave_speed_and_period as wavelength_law

# Description.
## Speed of light in air is 299910 km/s. Frequency of radio Europa+ is 101.6 MHz and therefore has oscillation period of 9,84252e-9 seconds. 
## According to online calculator https://vpayaem.ru/inf_wave1.html wavelength should be 2.95m.

@fixture
def test_args():
    v1 = Quantity(299910000 * units.meter / units.second)    
    period1 = Quantity(9.84252e-9 * units.second)
    Args = namedtuple("Args", ["v1", "period1"])
    return Args(v1=v1, period1=period1)

def test_basic_wavelength(test_args):
    result = wavelength_law.calculate_wavelength(test_args.v1, test_args.period1)    
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_wavelength = convert_to(result, units.meter).subs({units.meter: 1}).evalf(6)    
    assert result_wavelength == approx(2.95, 0.001)

def test_bad_velocity(test_args):
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wavelength_law.calculate_wavelength(vb, test_args.period1)
    with raises(TypeError):
        wavelength_law.calculate_wavelength(100, test_args.period1)

def test_bad_period(test_args):
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wavelength_law.calculate_wavelength(test_args.v1, pb)
    with raises(TypeError):
        wavelength_law.calculate_wavelength(test_args.v1, 100)
