# Description
## For example we are calculatiing gravity force for objects with masses 3000 and 5000 kilos
## with 0.06 meters of distance betweem their centers of mass.
## according to matematika-club.ru gravity calculator, force should be equal to 0.27809583 Newtons.

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
    SI.set_quantity_scale_factor(m1, 3000 * units.kilogram)
    m2 = units.Quantity('m2')
    SI.set_quantity_dimension(m2, units.mass)
    SI.set_quantity_scale_factor(m2, 5000 * units.kilogram)
    R = units.Quantity('R')
    SI.set_quantity_dimension(R, units.length)
    SI.set_quantity_scale_factor(R, 0.06 * units.meter)

    Args = namedtuple('Args', ['m1', 'm2', 'R'])
    return Args(m1 = m1, m2 = m2, R = R)

def test_basic_force(test_args):
    result = gravity_law.calculate_force(test_args.m1, test_args.m2, test_args.R)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)

    result_force = convert_to(result, units.newton).subs(units.newton, 1).evalf(7)
    assert result_force == approx(0.27809583, 0.0000001)


def test_bad_mass(test_args):
    mb = units.Quantity('mb')
    SI.set_quantity_dimension(mb, units.length)
    SI.set_quantity_scale_factor(mb, 1 * units.meter)

    with raises(errors.UnitsError):
        gravity_law.calculate_force(mb, test_args.m2, test_args.R)

    with raises(TypeError):
        gravity_law.calculate_force(100, test_args.m2, test_args.R)

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
