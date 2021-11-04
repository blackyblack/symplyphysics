from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.chemistry import avogadro_number_from_mole_count

@fixture
def test_args():
    mole_count = units.Quantity('mole_count')
    SI.set_quantity_dimension(mole_count, units.amount_of_substance)
    SI.set_quantity_scale_factor(mole_count, 5 * units.mole)

    Args = namedtuple('Args', ['M'])
    return Args(M = mole_count)

def test_basic_particles_count(test_args):
    result = avogadro_number_from_mole_count.calculate_particles_count(test_args.M)
    assert isinstance(result, int)
    assert result == approx(3.011E+24, 0.01)

def test_bad_mole_count():
    Mb = units.Quantity('Mb')
    SI.set_quantity_dimension(Mb, units.length)
    SI.set_quantity_scale_factor(Mb, 3 * units.meter)

    with raises(errors.UnitsError):
        avogadro_number_from_mole_count.calculate_particles_count(Mb)

    with raises(TypeError):
        avogadro_number_from_mole_count.calculate_particles_count(100)
