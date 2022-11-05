from collections import namedtuple
from pytest import approx, fixture

from math import pi
from sympy import sin, cos, sqrt
from sympy.vector import CoordSys3D, VectorZero
from symplyphysics import (
    units, convert_to, SI, expr_to_quantity
)
from symplyphysics.laws.fields import circulation_is_integral_of_curl_over_surface as circulation_def

@fixture
def test_args():
    C = CoordSys3D('C')
    force_unit = units.Quantity('force_unit')
    SI.set_quantity_dimension(force_unit, units.force)
    SI.set_quantity_scale_factor(force_unit, 1 * units.newton)
    radius_unit = units.Quantity('radius_unit')
    SI.set_quantity_dimension(radius_unit, units.length)
    SI.set_quantity_scale_factor(radius_unit, 1 * units.meter)
    Args = namedtuple('Args', ['C', 'force_unit', 'radius_unit'])
    return Args(C=C, force_unit=force_unit, radius_unit=radius_unit)

def test_basic_circulation(test_args):
    field = test_args.C.y * test_args.C.i + (test_args.C.z + test_args.C.x) * test_args.C.k
    surface = circulation_def.parameter1 * cos(circulation_def.parameter2) * test_args.C.i +  circulation_def.parameter1 * sin(circulation_def.parameter2) * test_args.C.j
    result_expr = circulation_def.calculate_circulation(field, surface, 0, 1, 0, pi / 2)
    assert result_expr.evalf(4) == approx(-pi / 4, 0.001)

def test_two_parameters_circulation(test_args):
    field = test_args.C.y * test_args.C.i + -1 * test_args.C.x * test_args.C.j
    # circle function is: x**2 + y**2 = 9
    # from circulation_is_integral_along_curve_test we got circulation -18 * pi
    # let's check with cone surface
    cone = 3 * circulation_def.parameter1 * cos(circulation_def.parameter2) * test_args.C.i + 3 * circulation_def.parameter1 * sin(circulation_def.parameter2) * test_args.C.j + circulation_def.parameter1 * test_args.C.k
    result_expr = circulation_def.calculate_circulation(field, cone, 0, 1, 0, 2 * pi)
    assert result_expr.evalf(4) == approx(-18 * pi, 0.001)

def test_gravitational_field_is_conservative(test_args):
    distance = sqrt(test_args.C.x**2 + test_args.C.y**2 + test_args.C.z**2)
    field = -test_args.C.x / (distance**2 * distance) * test_args.C.i - test_args.C.y / (distance**2 * distance) * test_args.C.j - test_args.C.z / (distance**2 * distance) * test_args.C.k
    field_rotor_applied = circulation_def.field_rotor_definition.rhs.subs(circulation_def.field, field).doit()
    assert field_rotor_applied == VectorZero.zero

def test_force_field_circulation(test_args):
    distance = sqrt(test_args.C.x**2 + test_args.C.y**2 + test_args.C.z**2)
    # we use gravitational force as reference
    # G = force * length**2 / mass**2
    #TODO: explain E
    #TODO: explain distance**3
    # E = G * m * M / r**2 = force * length**2 / r**2
    field_unit = (test_args.force_unit * test_args.radius_unit**2) / (distance**2 * distance)
    field = -test_args.C.y * field_unit * test_args.C.i + test_args.C.x * field_unit * test_args.C.j + test_args.C.z * field_unit * test_args.C.k
    surface = circulation_def.parameter1 * cos(circulation_def.parameter2) * test_args.C.i +  circulation_def.parameter1 * sin(circulation_def.parameter2) * test_args.C.j
    result_expr = circulation_def.calculate_circulation(field, surface, 1 * test_args.radius_unit, 2 * test_args.radius_unit, 0, pi / 2)
    result = expr_to_quantity(result_expr, 'force_work')
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_work = convert_to(result, units.joule).subs({units.joule: 1}).evalf(4)
    assert result_work == approx(-pi / 4, 0.001)
