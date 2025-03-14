from collections import namedtuple
from typing import Sequence, Any
from pytest import fixture, raises
from sympy import Expr, atan, cos, sin, sqrt, symbols
from sympy.vector import express
from symplyphysics.core.test_decorators import unsupported_usage
from symplyphysics.core.points.cylinder_point import CylinderPoint
from symplyphysics.core.points.cartesian_point import CartesianPoint
from symplyphysics.core.points.point import Point
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem, coordinates_rotate, coordinates_transform
from symplyphysics.core.fields.vector_field import VectorField


def _assert_point(field_: VectorField, point_: Point, expected_: Sequence[Any]) -> None:
    value = field_(point_)
    for idx, v in enumerate(value.components):
        assert v == expected_[idx]


Args = namedtuple("Args", ["C"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C = CoordinateSystem()
    return Args(C=C)


# Test VectorField constructor


def test_basic_field() -> None:

    def field_function(p: CartesianPoint) -> Sequence[Expr]:
        return [p.y, p.x]

    field = VectorField(field_function)
    assert callable(field.field_function)
    field_point = CartesianPoint(1, 2, 3)
    _assert_point(field, field_point, [2, 1])
    assert len(field.basis) == 3


def test_empty_field() -> None:
    field = VectorField([])
    assert not callable(field.field_function)
    field_point = Point(1, 2, 3)
    _assert_point(field, field_point, [])


def test_4d_field() -> None:
    field = VectorField(lambda p: [p.x, p.y, p.z, p.x])
    field_point = CartesianPoint(1, 2, 3)
    _assert_point(field, field_point, [1, 2, 3, 1])


def test_4d_point_field() -> None:
    field = VectorField(lambda p: [p.x, p.y, p.z, p.coordinate(3)])
    field_point = CartesianPoint(1, 2, 3, 4)
    _assert_point(field, field_point, [1, 2, 3, 4])


@unsupported_usage
def test_invalid_lambda_field() -> None:
    field = VectorField(lambda p: [p.y + "string", p.x])
    assert callable(field.field_function)
    field_point = CartesianPoint(1, 2, 3)
    # cannot add integer and string in field lambda
    with raises(TypeError):
        field(field_point)


def test_coord_system_field(test_args: Args) -> None:
    field = VectorField(lambda p: [p.y * p.z, 0, 0], test_args.C)
    field_point = CartesianPoint(1, 2, 3)
    _assert_point(field, field_point, [6, 0, 0])
    assert field.basis == [
        test_args.C.coord_system.x, test_args.C.coord_system.y, test_args.C.coord_system.z
    ]
    assert field.coordinate_system == test_args.C


# Test VectorField.from_sympy_vector()


def test_basic_vector_to_field_conversion(test_args: Args) -> None:
    field = VectorField.from_sympy_vector(
        test_args.C.coord_system.x * test_args.C.coord_system.i +
        test_args.C.coord_system.y * test_args.C.coord_system.j, test_args.C)
    assert callable(field.field_function)
    field_point = Point(1, 2, 3)
    _assert_point(field, field_point, [1, 2, 0])
    assert field.basis == [
        test_args.C.coord_system.x, test_args.C.coord_system.y, test_args.C.coord_system.z
    ]
    assert field.coordinate_system == test_args.C


def test_skip_dimension_vector_to_field_conversion(test_args: Args) -> None:
    field = VectorField.from_sympy_vector(
        1 * test_args.C.coord_system.i + 2 * test_args.C.coord_system.k, test_args.C)
    field_point = Point(1, 2, 3)
    _assert_point(field, field_point, [1, 0, 2])


# Base scalar is not a valid SymPy vector
def test_only_scalar_to_field_conversion(test_args: Args) -> None:
    with raises(AttributeError):
        VectorField.from_sympy_vector(test_args.C.coord_system.x, test_args.C)


# different coordinate systems in parameters are not supported
def test_different_coord_systems_vector_to_field_conversion(test_args: Args) -> None:
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    _, theta, _ = C1.coord_system.base_scalars()
    i, j, _ = C1.coord_system.base_vectors()
    with raises(TypeError):
        VectorField.from_sympy_vector(
            test_args.C.coord_system.x * test_args.C.coord_system.i + 2 * theta * j, test_args.C)
    with raises(TypeError):
        VectorField.from_sympy_vector(test_args.C.coord_system.x * i, test_args.C)


def test_custom_names_vector_to_field_conversion() -> None:
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    r, theta, z = C1.coord_system.base_scalars()
    i, j, _ = C1.coord_system.base_vectors()
    field = VectorField.from_sympy_vector(r * i + 2 * theta * j, C1)
    assert callable(field.field_function)
    field_point = Point(1, 2, 3)
    _assert_point(field, field_point, [1, 4, 0])
    assert field.basis == [r, theta, z]
    assert field.coordinate_system == C1


def test_rotate_coordinates_vector_to_field_conversion(test_args: Args) -> None:
    sympy_vector_field = test_args.C.coord_system.x * test_args.C.coord_system.i + test_args.C.coord_system.y * test_args.C.coord_system.j
    field = VectorField.from_sympy_vector(sympy_vector_field, test_args.C)
    assert callable(field.field_function)
    field_point = Point(1, 2, 3)
    _assert_point(field, field_point, [1, 2, 0])
    theta = symbols("theta")
    B = coordinates_rotate(test_args.C, theta, test_args.C.coord_system.k)
    x, y, z = B.coord_system.base_scalars()
    transformed_vector = express(sympy_vector_field, B.coord_system, variables=True)
    result_transformed_field = VectorField.from_sympy_vector(transformed_vector, B)
    _assert_point(result_transformed_field, field_point,
        [(sin(theta) + 2 * cos(theta)) * sin(theta) + (cos(theta) - 2 * sin(theta)) * cos(theta),
        (sin(theta) + 2 * cos(theta)) * cos(theta) - (cos(theta) - 2 * sin(theta)) * sin(theta), 0])
    assert result_transformed_field.basis == [x, y, z]
    assert result_transformed_field.coordinate_system == B


# when we express SymPy vector field to another coordinate system and base scalars (C.x, C.y) are
# not transformed, we get multiple coordinate systems, which is unsupported by sympy_vector_to_field
def test_rotate_coordinates_without_variables_vector_to_field_conversion(test_args: Args) -> None:
    theta = symbols("theta")
    B = coordinates_rotate(test_args.C, theta, test_args.C.coord_system.k)
    transformed_field = express(
        test_args.C.coord_system.x * test_args.C.coord_system.i +
        test_args.C.coord_system.y * test_args.C.coord_system.j, B.coord_system)
    with raises(TypeError):
        VectorField.from_sympy_vector(transformed_field, B)


# Test field.apply_to_basis()


def test_basic_apply_to_basis(test_args: Args) -> None:
    field = VectorField(lambda p: [p.y, p.x, 0], test_args.C)
    field_space = field.apply_to_basis()
    assert field_space.coordinate_system == test_args.C
    assert field_space.components == [test_args.C.coord_system.y, test_args.C.coord_system.x, 0]


def test_custom_names_apply_to_basis() -> None:

    def field_function(p: CylinderPoint) -> Sequence[Expr]:
        return [p.radius, p.azimuthal_angle, p.height]

    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    r, theta, z = C1.coord_system.base_scalars()
    field = VectorField(field_function, C1)
    field_space = field.apply_to_basis()
    assert field_space.coordinate_system == C1
    assert field_space.components == [r, theta, z]


# Test sympy_vector_from_field()


def test_basic_sympy_vector_from_field(test_args: Args) -> None:
    field = VectorField(lambda p: [p.y, p.x, 0], test_args.C)
    vector = field.to_sympy_vector()
    assert vector == test_args.C.coord_system.y * test_args.C.coord_system.i + test_args.C.coord_system.x * test_args.C.coord_system.j


def test_plain_sympy_vector_from_field(test_args: Args) -> None:
    field = VectorField([1, 2, 0], test_args.C)
    vector = field.to_sympy_vector()
    assert vector == 1 * test_args.C.coord_system.i + 2 * test_args.C.coord_system.j


# Test field.apply()


# Result is a function that returns a vector at any point of the trajectory.
# Result is stored in array, where first component of the array is magnitude of the resulting
# vector along X-axis (also called i-vector), second component is magnitude of the resulting
# vector along Y-axis (also called j-vector).
# Input field has X and Y swapped, so as we are moving along X-axis of the trajectory,
# resulting Y component of the vector grows.
def test_basic_field_apply(test_args: Args) -> None:
    field = VectorField(lambda p: [p.y, p.x], test_args.C)
    # represents surface
    trajectory = [test_args.C.coord_system.x, test_args.C.coord_system.y]
    trajectory_vectors = field.apply(trajectory)
    assert trajectory_vectors.coordinate_system == test_args.C
    assert trajectory_vectors.components == [test_args.C.coord_system.y, test_args.C.coord_system.x]


def test_parametrized_field_apply(test_args: Args) -> None:
    field = VectorField.from_sympy_vector(
        -test_args.C.coord_system.y * test_args.C.coord_system.i +
        test_args.C.coord_system.x * test_args.C.coord_system.j, test_args.C)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_vectors = field.apply(trajectory)
    assert trajectory_vectors.coordinate_system == test_args.C
    assert trajectory_vectors.components == [-1 * parameter, parameter, 0]


def test_sympy_field_apply(test_args: Args) -> None:
    field = VectorField.from_sympy_vector(
        -test_args.C.coord_system.y * test_args.C.coord_system.i +
        test_args.C.coord_system.x * test_args.C.coord_system.j, test_args.C)
    trajectory = [test_args.C.coord_system.x, test_args.C.coord_system.y]
    trajectory_vectors = field.apply(trajectory)
    assert trajectory_vectors.coordinate_system == test_args.C
    assert trajectory_vectors.components == [
        -test_args.C.coord_system.y, test_args.C.coord_system.x, 0
    ]


# VectorField is not rebased automatically and should be rebased to the same coordinate
# system as in trajectory with 'field_rebase'.
def test_different_coord_systems_field_apply(test_args: Args) -> None:
    result_field = VectorField(lambda p: [p.y * p.x, 0, 0], test_args.C)
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    r, theta, _ = C1.coord_system.base_scalars()
    trajectory = [r, theta]
    trajectory_vectors = result_field.apply(trajectory)
    assert trajectory_vectors.components == [r * theta, 0, 0]


# Test non-cartesian coordinate systems


def test_cylindrical_field_create(test_args: Args) -> None:
    field = VectorField(lambda p: [p.x, p.y, 0], test_args.C)
    point = [1, 2]
    point_vector = field.apply(point)
    assert point_vector.components == [1, 2, 0]

    B = coordinates_transform(test_args.C, CoordinateSystem.System.CYLINDRICAL)
    field_rebased = VectorField(lambda p: [p.r, p.theta, 0], B)
    assert field_rebased.coordinate_system == B

    # point should have r = sqrt(5) in polar coordinates
    # theta angle is atan(2/1)
    point_polar = [sqrt(5), atan(2)]
    point_polar_vector = field_rebased.apply(point_polar)
    assert point_polar_vector.components == [sqrt(5), atan(2), 0]

    # now rebase polar vector back to cartesian coordinates and confirm
    # it is the same as original vector
    vector_rebased = point_polar_vector.rebase(test_args.C)
    assert vector_rebased.components == [1, 2, 0]
