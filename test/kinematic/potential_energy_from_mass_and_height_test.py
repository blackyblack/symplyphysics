from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.kinematic import potential_energy_from_mass_and_height as potential_energy

@fixture
def test_args():
    m = units.Quantity('m')
    SI.set_quantity_dimension(m, units.mass)
    SI.set_quantity_scale_factor(m, 1 * units.kilogram)
    g = units.Quantity('g')
    SI.set_quantity_dimension(g, units.acceleration)
    SI.set_quantity_scale_factor(g, 9.82 * units.meter / units.second**2)
    h = units.Quantity('h')
    SI.set_quantity_dimension(h, units.length)
    SI.set_quantity_scale_factor(h, 1 * units.meter)
    Args = namedtuple('Args', ['m', 'g', 'h'])
    return Args(m=m, g=g, h=h)

def test_basic_energy(test_args):
    result = potential_energy.calculate_potential_energy(test_args.m, test_args.g, test_args.h)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).subs(units.joule, 1).evalf(3)
    assert result_energy == approx(9.82, 0.01)

def test_bad_body_mass(test_args):
    bm = units.Quantity('bm')
    SI.set_quantity_dimension(bm, units.acceleration)
    SI.set_quantity_scale_factor(bm, 1 * units.meter)
    with raises(errors.UnitsError):
        potential_energy.calculate_potential_energy(bm, test_args.g, test_args.h)
    with raises(TypeError):
        potential_energy.calculate_potential_energy(100, test_args.g, test_args.h)

def test_bad_gravitation_acceleration(test_args):
    bg = units.Quantity('bg')
    SI.set_quantity_dimension(bg, units.length)
    SI.set_quantity_scale_factor(bg, 1 * units.kilogram)
    with raises(errors.UnitsError):
        potential_energy.calculate_potential_energy(test_args.m, bg, test_args.h)
    with raises(TypeError):
        potential_energy.calculate_potential_energy(test_args.m, 100, test_args.h)

def test_bad_height(test_args):
    bh = units.Quantity('bh')
    SI.set_quantity_dimension(bh, units.acceleration)
    SI.set_quantity_scale_factor(bh, 1 * units.meter / units.second**2)
    with raises(errors.UnitsError):
        potential_energy.calculate_potential_energy(test_args.m, test_args.g, bh)
    with raises(TypeError):
        potential_energy.calculate_potential_energy(test_args.m, test_args.g, 100)