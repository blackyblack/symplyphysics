from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.dynamics import acceleration_from_force as newton_second_law

@fixture
def test_args():
    m = units.Quantity('m')
    SI.set_quantity_dimension(m, units.mass)
    SI.set_quantity_scale_factor(m, 1 * units.kilogram)
    a = units.Quantity('a')
    SI.set_quantity_dimension(a, units.acceleration)
    SI.set_quantity_scale_factor(a, 3 * units.meter / units.second**2)

    Args = namedtuple('Args', ['m', 'a'])
    return Args(m = m, a = a)

def test_basic_force(test_args):
    result = newton_second_law.calculate_force(test_args.m, test_args.a)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)

    result_force = convert_to(result, units.newton).subs(units.newton, 1).evalf(2)
    assert result_force == approx(3.0, 0.01)

def test_bad_mass(test_args):
    mb = units.Quantity('mb')
    SI.set_quantity_dimension(mb, units.length)
    SI.set_quantity_scale_factor(mb, 1 * units.meter)

    with raises(errors.UnitsError):
        newton_second_law.calculate_force(mb, test_args.a)

    with raises(TypeError):
        newton_second_law.calculate_force(100, test_args.a)

def test_bad_acceleration(test_args):
    ab = units.Quantity('ab')
    SI.set_quantity_dimension(ab, units.length)
    SI.set_quantity_scale_factor(ab, 3 * units.meter)

    with raises(errors.UnitsError):
        newton_second_law.calculate_force(test_args.m, ab)

    with raises(TypeError):
        newton_second_law.calculate_force(test_args.m, 100)
