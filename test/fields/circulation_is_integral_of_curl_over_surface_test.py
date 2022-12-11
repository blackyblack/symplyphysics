from collections import namedtuple
from pytest import approx, fixture

from math import pi
from sympy import Expr, sin, cos, sqrt
from symplyphysics import (
    units, SI, expr_to_quantity, convert_to,
    array_to_sympy_vector, CoordSys3D, VectorZero, FieldPoint, VectorField, apply_field
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
    field = VectorField(lambda point: point.y, 0, lambda point: point.x + point.z)
    surface = [circulation_def.parameter1 * cos(circulation_def.parameter2), circulation_def.parameter1 * sin(circulation_def.parameter2)]
    result_expr = circulation_def.calculate_circulation(test_args.C, field, surface, 0, 1, 0, pi / 2)
    assert result_expr.evalf(4) == approx(-pi / 4, 0.001)

def test_two_parameters_circulation(test_args):
    field = VectorField(lambda point: point.y, lambda point: -point.x)
    # circle function is: x**2 + y**2 = 9
    # from circulation_is_integral_along_curve_test we got circulation -18 * pi
    # let's check with cone surface
    cone = [3 * circulation_def.parameter1 * cos(circulation_def.parameter2), 3 * circulation_def.parameter1 * sin(circulation_def.parameter2), circulation_def.parameter1]
    result_expr = circulation_def.calculate_circulation(test_args.C, field, cone, 0, 1, 0, 2 * pi)
    assert result_expr.evalf(4) == approx(-18 * pi, 0.001)

def _distance(point: FieldPoint) -> Expr:
    return sqrt(point.x**2 + point.y**2 + point.z**2)

def test_gravitational_field_is_conservative(test_args):
    field = VectorField(lambda point: -point.x / _distance(point)**3, lambda point: -point.y / _distance(point)**3, lambda point: -point.z / _distance(point)**3)
    unit_coordinates = test_args.C.base_scalars()
    unit_trajectory = [unit_coordinates[0], unit_coordinates[1], unit_coordinates[2]]
    field_unit_app = apply_field(test_args.C, field, unit_trajectory)
    field_unit_sympy_vector = array_to_sympy_vector(test_args.C, field_unit_app)
    field_rotor_applied = circulation_def.field_rotor_definition.rhs.subs(field, field_unit_sympy_vector).doit()
    assert field_rotor_applied == VectorZero.zero

def test_force_field_circulation(test_args):
    # we use lorentz force in magnetic field as reference
    # B = mass / (current * time**2) = mass / (charge * time)
    # Lorentz force is: F = q * v * B = charge * (length / time) * B = force
    field = VectorField(lambda point: -point.y * test_args.force_unit / _distance(point), lambda point: point.x * test_args.force_unit / _distance(point))
    surface = [circulation_def.parameter1 * cos(circulation_def.parameter2), circulation_def.parameter1 * sin(circulation_def.parameter2)]
    result_expr = circulation_def.calculate_circulation(test_args.C, field, surface, 1 * test_args.radius_unit, 2 * test_args.radius_unit, 0, pi / 2)
    result = expr_to_quantity(result_expr, 'force_work')
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_work = convert_to(result, units.joule).subs({units.joule: 1}).evalf(2)
    assert result_work > 0
