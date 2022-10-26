from collections import namedtuple
from pytest import approx, fixture, raises

from math import pi
from sympy import sin, cos
from sympy.vector import CoordSys3D
from symplyphysics import (
    units, convert_to, SI, errors, expr_to_quantity
)
from symplyphysics.definitions import circulation_is_integral_along_curve as circulation_def

@fixture
def test_args():
    force_unit = units.Quantity('force_unit')
    SI.set_quantity_dimension(force_unit, units.force)
    SI.set_quantity_scale_factor(force_unit, 1 * units.newton)
    radius_unit = units.Quantity('radius_unit')
    SI.set_quantity_dimension(radius_unit, units.length)
    SI.set_quantity_scale_factor(radius_unit, 1 * units.meter)

    C = CoordSys3D('C')
    # trajectory is ellipse: x^2 / 4 + y^2 / 9 = 1
    # parametrize by t (use polar coordinates): x(t) = 2 * cos(t), y(t) = 3 * sin(t)
    trajectory_x = 2 * cos(circulation_def.parameter)
    trajectory_y = 3 * sin(circulation_def.parameter)
    trajectory = trajectory_x * radius_unit * C.i + trajectory_y * radius_unit * C.j
    # field is a field of forces
    # field is (-y, x): -1 * C.y * C.i + C.x * C.j
    # parametrize by t
    field = -1 * trajectory_y * force_unit * C.i + trajectory_x * force_unit * C.j

    # we parametrized with polar coordinates so integrate full circle [0, 2 * pi]
    position_from = 0
    position_to = 2 * pi

    Args = namedtuple('Args', ['field', 'trajectory', 'p0', 'p1'])
    return Args(field = field, trajectory = trajectory, p0 = position_from, p1 = position_to)

def test_force_circulation(test_args):
    result_expr = circulation_def.calculate_circulation(test_args.field, test_args.trajectory, test_args.p0, test_args.p1)
    result = expr_to_quantity(result_expr, 'force_work')
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_work = convert_to(result, units.joule).subs({units.joule: 1}).evalf(4)
    assert result_work == approx(12 * pi, 0.001)

'''
def test_bad_mass(test_args):
    mb = units.Quantity('mb')
    SI.set_quantity_dimension(mb, units.length)
    SI.set_quantity_scale_factor(mb, 1 * units.meter)

    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(mb, test_args.v)

    with raises(TypeError):
        momentum_def.calculate_momentum(100, test_args.v)

def test_bad_velocity(test_args):
    vb = units.Quantity('vb')
    SI.set_quantity_dimension(vb, units.length)
    SI.set_quantity_scale_factor(vb, 1 * units.meter)

    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(test_args.m, vb)

    with raises(TypeError):
        momentum_def.calculate_momentum(test_args.m, 100)
'''
