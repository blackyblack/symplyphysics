from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.definitions import refractive_index_is_wave_speeds_ratio as refractive_index_definition

# Description.
## Known relative refractive index of air-water is 1.298. Propagation speed of visible light in air is 299910 km/s.
## Propagation speed in water is 231000 km/s.

@fixture
def test_args():
    v1 = units.Quantity('v1')
    SI.set_quantity_dimension(v1, units.velocity)
    SI.set_quantity_scale_factor(v1, 299910 * units.kilometer / units.second)
        
    v2 = units.Quantity('v2')
    SI.set_quantity_dimension(v2, units.velocity)
    SI.set_quantity_scale_factor(v2, 231000 * units.kilometer / units.second)

    Args = namedtuple('Args', ['v1', 'v2'])
    return Args(v1=v1, v2=v2)

def test_basic_refraction_factor(test_args):
    result = refractive_index_definition.calculate_refractive_index(test_args.v1, test_args.v2)    
    assert result == approx(1.298, 0.001)

def test_bad_velocity(test_args):
    vb = units.Quantity('vb')
    SI.set_quantity_dimension(vb, units.length)
    SI.set_quantity_scale_factor(vb, 1 * units.meter)
    with raises(errors.UnitsError):
        refractive_index_definition.calculate_refractive_index(vb, test_args.v2)

    with raises(TypeError):
        refractive_index_definition.calculate_refractive_index(100, test_args.v2)

    with raises(errors.UnitsError):
        refractive_index_definition.calculate_refractive_index(test_args.v1, vb)

    with raises(TypeError):
        refractive_index_definition.calculate_refractive_index(test_args.v1, 100)

