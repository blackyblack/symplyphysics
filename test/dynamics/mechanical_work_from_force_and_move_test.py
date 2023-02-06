from collections import namedtuple
from pytest import approx, fixture
from symplyphysics import (
    units, convert_to, SI
)
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type
from symplyphysics.laws.dynamics import mechanical_work_from_force_and_move as work_law

@fixture
def test_args():
    F = units.Quantity('F')
    SI.set_quantity_dimension(F, units.force)
    SI.set_quantity_scale_factor(F, 100 * units.newton)
    S = units.Quantity('S')
    SI.set_quantity_dimension(S, units.length)
    SI.set_quantity_scale_factor(S, 3 * units.meter)
    # Force is applied at 60 degrees angle but body moves horisontally
    angle_between_force_and_axis = units.Quantity('angle_between_force_and_axis')
    SI.set_quantity_dimension(angle_between_force_and_axis, angle_type)
    SI.set_quantity_scale_factor(angle_between_force_and_axis, 60 * units.degree)
    angle_between_movement_and_axis = units.Quantity('angle_between_movement_and_axis')
    SI.set_quantity_dimension(angle_between_movement_and_axis, angle_type)
    SI.set_quantity_scale_factor(angle_between_movement_and_axis, 0 * units.degree)

    Args = namedtuple('Args', ['F', 'S', 'Fa', 'Sa'])
    return Args(F = F, S = S, Fa = angle_between_force_and_axis, Sa = angle_between_movement_and_axis)

def test_basic_work(test_args):
    result = work_law.calculate_work(test_args.F, test_args.S, test_args.Fa, test_args.Sa)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)

    result_work = convert_to(result, units.joule).subs(units.joule, 1).evalf(4)
    assert result_work == approx(150, 0.01)

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