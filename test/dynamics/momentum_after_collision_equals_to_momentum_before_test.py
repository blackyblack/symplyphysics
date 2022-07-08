from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)

from symplyphysics.laws.dynamics import momentum_after_collision_equals_to_momentum_before as momentum_law

@fixture
def test_args():
    P_before = units.Quantity('P_before')
    SI.set_quantity_dimension(P_before, units.momentum)
    SI.set_quantity_scale_factor(P_before, 5 * units.kilogram * units.meter / units.second)   

    Args = namedtuple('Args', ['P_before'])
    return Args(P_before = P_before)

def test_basic_functionality(test_args):
    result = momentum_law.calculate_momentum_after(test_args.P_before)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.momentum)

    result_ = convert_to(result, units.momentum).subs(units.kilogram * units.meter / units.second, 1).evalf(2)
    assert result_ == approx(5.0, 0.01)

"""
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
"""