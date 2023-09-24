from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import S, Expr, sin, cos, sqrt, pi
from sympy.vector import VectorZero
from symplyphysics import (
    Symbol,
    units,
    Quantity,
    CoordinateSystem,
    SI,
    convert_to,
)
from symplyphysics.core.fields.field_point import FieldPoint
from symplyphysics.core.fields.operators import curl_operator
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.laws.fields import circulation_is_integral_of_curl_over_surface as circulation_def


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    force_unit = Quantity(1 * units.newton)
    radius_unit = Quantity(1 * units.meter)
    Args = namedtuple("Args", ["C", "force_unit", "radius_unit"])
    return Args(C=C, force_unit=force_unit, radius_unit=radius_unit)


def test_basic_circulation(test_args):
    field = VectorField(lambda point: [point.y, 0, point.x + point.z], test_args.C)
    surface = [
        circulation_def.parameter1 * cos(circulation_def.parameter2),
        circulation_def.parameter1 * sin(circulation_def.parameter2)
    ]
    result = circulation_def.calculate_circulation(field, surface,
        [circulation_def.parameter1, circulation_def.parameter2], (0, 1), (0, pi / 2))
    assert convert_to(result, S.One).evalf(4) == approx((-pi / 4).evalf(4), 0.001)


def test_two_parameters_circulation(test_args):
    field = VectorField(lambda point: [point.y, -point.x, 0], test_args.C)
    # circle function is: x**2 + y**2 = 9
    # from circulation_is_integral_along_curve_test we got circulation -18 * pi
    # let's check with cone surface
    cone = [
        3 * circulation_def.parameter1 * cos(circulation_def.parameter2),
        3 * circulation_def.parameter1 * sin(circulation_def.parameter2), circulation_def.parameter1
    ]
    result = circulation_def.calculate_circulation(field, cone,
        [circulation_def.parameter1, circulation_def.parameter2], (0, 1), (0, 2 * pi))
    assert convert_to(result, S.One).evalf(4) == approx((-18 * pi).evalf(4), 0.001)


def _distance(point: FieldPoint) -> Expr:
    return sqrt(point.x**2 + point.y**2 + point.z**2)


def test_gravitational_field_is_conservative(test_args):
    # gravitational field also has a common multiplier of -G*m1*m2. It does not
    # affect conservative property of a field.
    field = VectorField(
        lambda point: [
        point.x / _distance(point)**3, point.y / _distance(point)**3, point.z / _distance(point)**3
        ], test_args.C)
    field_rotor = curl_operator(field)
    field_rotor_applied = field_rotor.apply_to_basis().to_sympy_vector()
    assert field_rotor_applied == VectorZero.zero


def test_force_field_circulation(test_args):
    # we use lorentz force in magnetic field as reference
    # B = mass / (current * time**2) = mass / (charge * time)
    # Lorentz force is: F = q * v * B = charge * (length / time) * B
    field = VectorField(
        lambda point: [
        -point.y * test_args.force_unit / _distance(point), point.x * test_args.force_unit /
        _distance(point), 0
        ], test_args.C)
    surface = [
        circulation_def.parameter1 * cos(circulation_def.parameter2),
        circulation_def.parameter1 * sin(circulation_def.parameter2)
    ]
    result = circulation_def.calculate_circulation(field, surface,
        [circulation_def.parameter1, circulation_def.parameter2],
        (1 * test_args.radius_unit, 2 * test_args.radius_unit), (0, pi / 2))
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_work = convert_to(result, units.joule).evalf(2)
    assert result_work > 0


def test_three_parameters_circulation(test_args):
    field = VectorField(lambda point: [point.y, -point.x, 0], test_args.C)
    parameter3 = Symbol("parameter3")
    surface = [
        circulation_def.parameter1 * cos(circulation_def.parameter2),
        circulation_def.parameter1 * sin(circulation_def.parameter2)
    ]
    with raises(ValueError):
        circulation_def.calculate_circulation(field, surface,
            [circulation_def.parameter1, circulation_def.parameter2, parameter3], (0, 1),
            (0, 2 * pi))
    with raises(ValueError):
        circulation_def.calculate_circulation(field, surface,
            [circulation_def.parameter1, circulation_def.parameter2, circulation_def.parameter1],
            (0, 1), (0, 2 * pi))
