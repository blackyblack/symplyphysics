from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.dynamics import buoyant_force_from_density_and_volume as archimedes_law

@fixture
def test_args():
    pf = units.Quantity('pf')
    SI.set_quantity_dimension(pf, units.mass / units.volume)
    # water density
    SI.set_quantity_scale_factor(pf, 1000 * units.kilogram / units.meter**3)
    V = units.Quantity('V')
    SI.set_quantity_dimension(V, units.volume)
    SI.set_quantity_scale_factor(V, 0.2 * units.meter**3)

    Args = namedtuple('Args', ['V', 'pf'])
    return Args(V = V, pf = pf)

def test_basic_force(test_args):
    result = archimedes_law.calculate_force_buoyant(test_args.pf, test_args.V)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)

    result_force = convert_to(result, units.newton).subs(units.newton, 1).evalf(4)
    assert result_force == approx(1961.3, 0.01)

def test_bad_density(test_args):
    pb = units.Quantity('pb')
    SI.set_quantity_dimension(pb, units.length)
    SI.set_quantity_scale_factor(pb, 1 * units.meter)

    with raises(errors.UnitsError):
        archimedes_law.calculate_force_buoyant(pb, test_args.V)

    with raises(TypeError):
        archimedes_law.calculate_force_buoyant(100, test_args.V)

def test_bad_volume(test_args):
    Vb = units.Quantity('Vb')
    SI.set_quantity_dimension(Vb, units.length)
    SI.set_quantity_scale_factor(Vb, 1 * units.meter)

    with raises(errors.UnitsError):
        archimedes_law.calculate_force_buoyant(test_args.pf, Vb)

    with raises(TypeError):
        archimedes_law.calculate_force_buoyant(test_args.pf, 100)
