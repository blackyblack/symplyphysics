from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.gravity import gravity_force_from_mass as gravity_law

@fixture
def test_args():
    m1 = units.Quantity('m1')
    SI.set_quantity_dimension(m1, units.mass)
    SI.set_quantity_scale_factor(m1, 20000 * units.kilogram)
    m2 = units.Quantity('m2')
    SI.set_quantity_dimension(m2, units.mass)
    SI.set_quantity_scale_factor(m2, 10000 * units.kilogram)
    R = units.Quantity('R')
    SI.set_quantity_dimension(R, units.length)
    SI.set_quantity_scale_factor(R, 1 * units.meter)

    Args = namedtuple('Args', ['m1', 'm2', 'R'])
    return Args(m1 = m1, m2 = m2, R = R)

def test_basic_force(test_args):
    result = gravity_law.calculate_force(test_args.m1, test_args.m2, test_args.R)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)

    result_force = convert_to(result, units.newton).subs(units.newton, 1).evalf(7)
    assert result_force == approx(0.0133486, 0.0000001)


def test_bad_mass_1(test_args):
    mb = units.Quantity('mb')
    SI.set_quantity_dimension(mb, units.length)
    SI.set_quantity_scale_factor(mb, 1 * units.meter)

    with raises(errors.UnitsError):
        gravity_law.calculate_force(mb, test_args.m2, test_args.R)

    with raises(TypeError):
        gravity_law.calculate_force(100, test_args.m2, test_args.R)


def test_bad_mass_2(test_args):
    mb = units.Quantity('mb')
    SI.set_quantity_dimension(mb, units.length)
    SI.set_quantity_scale_factor(mb, 1 * units.meter)

    with raises(errors.UnitsError):
        gravity_law.calculate_force(test_args.m1, mb, test_args.R)

    with raises(TypeError):
        gravity_law.calculate_force(test_args.m1, 100, test_args.R)        


def test_bad_distance(test_args):
    db = units.Quantity('db')
    SI.set_quantity_dimension(db, units.charge)
    SI.set_quantity_scale_factor(db, 1 * units.coulomb)

    with raises(errors.UnitsError):
        gravity_law.calculate_force(test_args.m1, test_args.m2, db)

    with raises(TypeError):
        gravity_law.calculate_force(test_args.m1, test_args.m2, 100)        
