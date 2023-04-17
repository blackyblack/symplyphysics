from collections import namedtuple
from pytest import approx, fixture

from sympy import Expr, sin, cos, sqrt, pi
from symplyphysics import (
    units, SI, expr_to_quantity, convert_to,
    VectorZero, FieldPoint, VectorField, CoordinateSystem, sympy_vector_from_vector
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.fields import circulation_is_integral_of_curl_over_surface as circulation_def

@fixture
def test_args():
    C = CoordinateSystem()
    force_unit = Quantity(1 * units.newton)
    radius_unit = Quantity(1 * units.meter)
    Args = namedtuple("Args", ["C", "force_unit", "radius_unit"])
    return Args(C=C, force_unit=force_unit, radius_unit=radius_unit)

def test_basic_circulation(test_args):
    field = VectorField(lambda point: point.y, lambda point: 0, lambda point: point.x + point.z, test_args.C)
    surface = [circulation_def.parameter1 * cos(circulation_def.parameter2), circulation_def.parameter1 * sin(circulation_def.parameter2)]
    result_expr = circulation_def.calculate_circulation(field, surface, 0, 1, 0, pi/2)
    assert result_expr.evalf(4) == approx((-pi / 4).evalf(4), 0.001)

def test_two_parameters_circulation(test_args):
    field = VectorField(lambda point: point.y, lambda point: -point.x, 0, test_args.C)
    # circle function is: x**2 + y**2 = 9
    # from circulation_is_integral_along_curve_test we got circulation -18 * pi
    # let's check with cone surface
    cone = [3 * circulation_def.parameter1 * cos(circulation_def.parameter2), 3 * circulation_def.parameter1 * sin(circulation_def.parameter2), circulation_def.parameter1]
    result_expr = circulation_def.calculate_circulation(field, cone, 0, 1, 0, 2 * pi)
    assert result_expr.evalf(4) == approx((-18 * pi).evalf(4), 0.001)

def _distance(point: FieldPoint) -> Expr:
    return sqrt(point.x**2 + point.y**2 + point.z**2)

def test_gravitational_field_is_conservative(test_args):
    field = VectorField(lambda point: -point.x / _distance(point)**3, lambda point: -point.y / _distance(point)**3, lambda point: -point.z / _distance(point)**3, test_args.C)
    field_space = field.apply_to_basis()
    field_space_sympy = sympy_vector_from_vector(field_space)
    field_rotor_applied = circulation_def.field_rotor_definition.rhs.subs(circulation_def.field, field_space_sympy).doit()
    assert field_rotor_applied == VectorZero.zero

def test_force_field_circulation(test_args):
    # we use lorentz force in magnetic field as reference
    # B = mass / (current * time**2) = mass / (charge * time)
    # Lorentz force is: F = q * v * B = charge * (length / time) * B = force
    field = VectorField(lambda point: -point.y * test_args.force_unit / _distance(point), lambda point: point.x * test_args.force_unit / _distance(point), 0, test_args.C)
    surface = [circulation_def.parameter1 * cos(circulation_def.parameter2), circulation_def.parameter1 * sin(circulation_def.parameter2)]
    result_expr = circulation_def.calculate_circulation(field, surface, 1 * test_args.radius_unit, 2 * test_args.radius_unit, 0, pi / 2)
    result = expr_to_quantity(result_expr)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_work = convert_to(result, units.joule).subs(units.joule, 1).evalf(2)
    assert result_work > 0
