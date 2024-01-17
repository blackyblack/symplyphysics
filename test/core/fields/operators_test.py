from collections import namedtuple
from typing import Sequence
from pytest import fixture
from sympy import Expr, cos, exp, sin, Symbol as SymSymbol, sqrt
from sympy.vector import VectorZero
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.fields.scalar_field import ScalarField
from symplyphysics.core.fields.operators import curl_operator, divergence_operator, gradient_operator
from symplyphysics.core.points.cartesian_point import CartesianPoint
from symplyphysics.core.points.cylinder_point import CylinderPoint
from symplyphysics.core.points.sphere_point import SpherePoint

# Tests are mostly based on vector calculus slides: https://www.slideshare.net/garghanish/coordinate-systems-and-vector-calculus


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    parameter1 = SymSymbol("parameter1")
    parameter2 = SymSymbol("parameter2")
    Args = namedtuple("Args", ["C", "parameter1", "parameter2"])
    return Args(C=C, parameter1=parameter1, parameter2=parameter2)


def test_basic_gradient(test_args):
    field = ScalarField(
        lambda point: point.x**2 + point.y**2 - point.z**2
    )
    result_vector = gradient_operator(field)

    x, y, z = field.coordinate_system.coord_system.base_scalars()
    correct_result = [2 * x, 2 * y, -2 * z]

    for result_component, correct_component in zip(result_vector.components, correct_result):
        assert expr_equals(result_component, correct_component)


def test_cylindrical_gradient(test_args):
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)

    def field_function(p: CylinderPoint) -> ScalarValue:
        return p.r ** 2 + p.r * p.z * sin(p.theta) - p.z ** 2

    cylindrical_field = ScalarField(field_function, C1)
    result_cylindrical_vector = gradient_operator(cylindrical_field)

    r, theta, z = cylindrical_field.coordinate_system.coord_system.base_scalars()
    correct_result = [
        2 * r + z * sin(theta), 
        z * cos(theta), 
        -2 * z + r * sin(theta),
    ]
    
    for result_component, correct_component in zip(result_cylindrical_vector.components, correct_result):
        assert expr_equals(result_component, correct_component)


def test_spherical_gradient(test_args):
    C1 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)

    def field_function(p: SpherePoint) -> ScalarValue:
        return - p.r**2 * cos(2 * p.phi + p.theta)

    spherical_field = ScalarField(field_function, C1)
    result_spherical_vector = gradient_operator(spherical_field)

    r, theta, phi = spherical_field.coordinate_system.coord_system.base_scalars()
    correct_result = [
        -2 * r * cos(2*phi + theta),
        r * sin(2*phi + theta) / sin(phi),
        2 * r * sin(2*phi + theta),
    ]

    for result_component, correct_component in zip(result_spherical_vector.components, correct_result):
        assert expr_equals(result_component, correct_component)


def test_basic_divergence(test_args):
    field = VectorField(
        lambda point: [exp(point.x) * cos(point.y),
        exp(point.x) * sin(point.y), point.z], test_args.C)
    result = divergence_operator(field)
    x = field.coordinate_system.coord_system.base_scalars()[0]
    y = field.coordinate_system.coord_system.base_scalars()[1]
    assert result == 2 * exp(x) * cos(y) + 1


def test_cylindrical_divergence(test_args):
    field = VectorField(lambda point: [point.x * 2 / 3, point.y * 2 / 3, point.z * 2 / 3],
        test_args.C)
    result = divergence_operator(field)
    assert result == 2
    # verify that in cylindrical coordinates result is same
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)

    def field_function(p: CylinderPoint) -> Sequence[ScalarValue]:
        return [p.radius * 2 / 3, 0, p.height * 2 / 3]

    cylindrical_field = VectorField(field_function, C1)
    result = divergence_operator(cylindrical_field)
    assert result == 2


def test_spherical_divergence():
    C1 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)

    def field_function(p: SpherePoint) -> Sequence[ScalarValue]:
        return [
            1 / p.radius**2 * cos(p.polar_angle),
            cos(p.polar_angle), p.radius * sin(p.polar_angle) * cos(p.azimuthal_angle)
        ]

    field = VectorField(field_function, C1)
    result = divergence_operator(field)
    theta = field.coordinate_system.coord_system.base_scalars()[1]
    phi = field.coordinate_system.coord_system.base_scalars()[2]
    assert expr_equals(result, 2 * cos(theta) * cos(phi))


def test_basic_curl(test_args):
    field = VectorField(lambda point: [point.x**2 * point.y * point.z, 0, point.x * point.z],
        test_args.C)
    result_field = curl_operator(field)
    x = field.coordinate_system.coord_system.base_scalars()[0]
    y = field.coordinate_system.coord_system.base_scalars()[1]
    z = field.coordinate_system.coord_system.base_scalars()[2]
    result_vector = result_field.apply_to_basis()
    assert expr_equals(result_vector.components[0], 0)
    assert expr_equals(result_vector.components[1], x**2 * y - z)
    assert expr_equals(result_vector.components[2], -x**2 * z)


def test_cylindrical_curl():
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)

    def field_function(p: CylinderPoint) -> Sequence[ScalarValue]:
        return [
            p.radius * sin(p.azimuthal_angle), p.radius**2 * p.height,
            p.height * cos(p.azimuthal_angle)
        ]

    field = VectorField(field_function, C1)
    result_field = curl_operator(field)
    r = field.coordinate_system.coord_system.base_scalars()[0]
    theta = field.coordinate_system.coord_system.base_scalars()[1]
    z = field.coordinate_system.coord_system.base_scalars()[2]
    result_vector = result_field.apply_to_basis()
    assert expr_equals(result_vector.components[0], -1 / r * (z * sin(theta) + r**3))
    assert expr_equals(result_vector.components[1], 0)
    assert expr_equals(result_vector.components[2], 3 * r * z - cos(theta))


def test_spherical_curl():
    C1 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)

    def field_function(p: SpherePoint) -> Sequence[ScalarValue]:
        return [
            1 / p.radius**2 * cos(p.polar_angle),
            cos(p.polar_angle), p.radius * sin(p.polar_angle) * cos(p.azimuthal_angle)
        ]

    field = VectorField(field_function, C1)
    result_field = curl_operator(field)
    r = field.coordinate_system.coord_system.base_scalars()[0]
    theta = field.coordinate_system.coord_system.base_scalars()[1]
    phi = field.coordinate_system.coord_system.base_scalars()[2]
    result_vector = result_field.apply_to_basis()
    assert expr_equals(result_vector.components[0], cos(2 * phi) / (r * sin(phi)) + sin(theta))
    assert expr_equals(result_vector.components[1], (2 * cos(theta) + 1 / r**3) * sin(phi))
    assert expr_equals(result_vector.components[2], -cos(phi) / r)


def _distance(point: CartesianPoint) -> Expr:
    return sqrt(point.x**2 + point.y**2 + point.z**2)


def test_gravitational_field_is_conservative(test_args):
    # gravitational field also has a common multiplier of -G*M. It does not
    # affect conservative property of a field.
    field = VectorField(
        lambda point: [
        point.x / _distance(point)**3, point.y / _distance(point)**3, point.z / _distance(point)**3
        ], test_args.C)
    field_rotor = curl_operator(field)
    field_rotor_applied = field_rotor.apply_to_basis().to_sympy_vector()
    assert field_rotor_applied == VectorZero.zero
