'''
def calculate_force(generator_mass_: Quantity, object_mass_: Quantity, r_distance_: Quantity) -> Quantity:
    result_force_expr = solve(law, gravity_force, dict=True)[0][gravity_force]
    result_expr = result_force_expr.subs({generator_mass: generator_mass_, object_mass: object_mass_, r_distance: r_distance_})
    return expr_to_quantity(result_expr, 'gravity_force')
'''

from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.gravity import gravity_force_from_mass as gravity_law

@fixture
def test_args():
    M = units.Quantity('M')
    SI.set_quantity_dimension(M, units.mass)
    SI.set_quantity_scale_factor(M, 1 * units.kilogram)
    m = units.Quantity('m')
    SI.set_quantity_dimension(m, units.mass)
    SI.set_quantity_scale_factor(m, 0.1 * units.kilogram)
    R = units.Quantity('R')
    SI.set_quantity_dimension(R, units.length)
    SI.set_quantity_scale_factor(R, 0.1 * units.kilogram)

    Args = namedtuple('Args', ['M', 'm', 'R'])
    return Args(M = M, m = m, R = R)

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
