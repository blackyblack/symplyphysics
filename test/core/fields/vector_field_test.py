from collections import namedtuple
from typing import Sequence
from test.test_decorators import unsupported_usage
from pytest import fixture, raises
from sympy import Expr, atan, cos, pi, sin, sqrt, symbols
from sympy.vector import express
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem, coordinates_rotate, coordinates_transform
from symplyphysics.core.vectors.vectors import vector_rebase
from symplyphysics.core.fields.field_point import FieldPoint
from symplyphysics.core.fields.vector_field import VectorField, field_from_sympy_vector, field_rebase, sympy_vector_from_field


def _assert_point(field_: VectorField, point_: FieldPoint, expected_: Sequence[Expr | float]):
    value = field_(point_)
    for idx, v in enumerate(value.components):
        assert v == expected_[idx]


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    Args = namedtuple("Args", ["C"])
    return Args(C=C)


# Test VectorField constructor


def test_basic_field(test_args):
    field = VectorField(test_args.C, lambda p: [p.y, p.x])
    assert callable(field.field_function)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [2, 1])
    assert len(field.basis) == 3


def test_empty_field(test_args):
    field = VectorField(test_args.C, [])
    assert not callable(field.field_function)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [])


def test_4d_field(test_args):
    field = VectorField(test_args.C, lambda p: [p.x, p.y, p.z, p.x])
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [1, 2, 3, 1])


def test_4d_point_field(test_args):
    field = VectorField(test_args.C, lambda p: [p.x, p.y, p.z, p.coordinate(3)])
    field_point = FieldPoint(1, 2, 3, 4)
    _assert_point(field, field_point, [1, 2, 3, 4])


@unsupported_usage
def test_wrong_type_lambda_field(test_args):
    field = VectorField(test_args.C, lambda p: ["string", p.x])
    assert callable(field.field_function)
    field_point = FieldPoint(1, 2, 3)
    # non expression lambda in a field is not processed and returns as is
    _assert_point(field, field_point, ["string", 1])


@unsupported_usage
def test_wrong_type_value_field(test_args):
    field = VectorField(test_args.C, ["string", 1])
    assert not callable(field.field_function)
    field_point = FieldPoint(1, 2, 3)
    # non expression in a field is not processed and returns as is
    _assert_point(field, field_point, ["string", 1])


@unsupported_usage
def test_invalid_lambda_field(test_args):
    field = VectorField(test_args.C, lambda p: [p.y + "string", p.x])
    assert callable(field.field_function)
    field_point = FieldPoint(1, 2, 3)
    # cannot add integer and string in field lambda
    with raises(TypeError):
        field(field_point)


@unsupported_usage
def test_effect_in_lambda_field(test_args):
    field = VectorField(test_args.C, lambda p: [f"{p.y}", p.x])
    assert callable(field.field_function)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, ["2", 1])


def test_coord_system_field(test_args):
    field = VectorField(test_args.C, lambda p: [p.y * p.z, 0, 0])
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [6, 0, 0])
    assert field.basis == [
        test_args.C.coord_system.x, test_args.C.coord_system.y, test_args.C.coord_system.z
    ]
    assert field.coordinate_system == test_args.C


# Test field_from_sympy_vector()


def test_basic_vector_to_field_conversion(test_args):
    field = field_from_sympy_vector(
        test_args.C.coord_system.x * test_args.C.coord_system.i +
        test_args.C.coord_system.y * test_args.C.coord_system.j, test_args.C)
    assert callable(field.field_function)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [1, 2, 0])
    assert field.basis == [
        test_args.C.coord_system.x, test_args.C.coord_system.y, test_args.C.coord_system.z
    ]
    assert field.coordinate_system == test_args.C


def test_skip_dimension_vector_to_field_conversion(test_args):
    field = field_from_sympy_vector(1 * test_args.C.coord_system.i + 2 * test_args.C.coord_system.k,
        test_args.C)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [1, 0, 2])


# Integer is not a valid SymPy vector
def test_only_integer_vector_to_field_conversion(test_args):
    with raises(AttributeError):
        field_from_sympy_vector(1, test_args.C)


# Base scalar is not a valid SymPy vector
def test_only_scalar_to_field_conversion(test_args):
    with raises(AttributeError):
        field_from_sympy_vector(test_args.C.coord_system.x, test_args.C)


# different coordinate systems in parameters are not supported
def test_different_coord_systems_vector_to_field_conversion(test_args):
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    with raises(TypeError):
        field_from_sympy_vector(
            test_args.C.coord_system.x * test_args.C.coord_system.i +
            2 * C1.coord_system.theta * C1.coord_system.j, test_args.C)
    with raises(TypeError):
        field_from_sympy_vector(test_args.C.coord_system.x * C1.coord_system.i, test_args.C)


def test_custom_names_vector_to_field_conversion():
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    field = field_from_sympy_vector(
        C1.coord_system.r * C1.coord_system.i + 2 * C1.coord_system.theta * C1.coord_system.j, C1)
    assert callable(field.field_function)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [1, 4, 0])
    assert field.basis == [C1.coord_system.r, C1.coord_system.theta, C1.coord_system.z]
    assert field.coordinate_system == C1


def test_rotate_coordinates_vector_to_field_conversion(test_args):
    sympy_vector_field = test_args.C.coord_system.x * test_args.C.coord_system.i + test_args.C.coord_system.y * test_args.C.coord_system.j
    field = field_from_sympy_vector(sympy_vector_field, test_args.C)
    assert callable(field.field_function)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [1, 2, 0])
    theta = symbols("theta")
    B = coordinates_rotate(test_args.C, theta, test_args.C.coord_system.k)
    transformed_vector = express(sympy_vector_field, B.coord_system, variables=True)
    result_transformed_field = field_from_sympy_vector(transformed_vector, B)
    _assert_point(result_transformed_field, field_point,
        [(sin(theta) + 2 * cos(theta)) * sin(theta) + (cos(theta) - 2 * sin(theta)) * cos(theta),
        (sin(theta) + 2 * cos(theta)) * cos(theta) - (cos(theta) - 2 * sin(theta)) * sin(theta), 0])
    assert result_transformed_field.basis == [B.coord_system.x, B.coord_system.y, B.coord_system.z]
    assert result_transformed_field.coordinate_system == B


# when we express SymPy vector field to another coordinate system and base scalars (C.x, C.y) are
# not transformed, we get multiple coordinate systems, which is unsupported by sympy_vector_to_field
def test_rotate_coordinates_without_variables_vector_to_field_conversion(test_args):
    theta = symbols("theta")
    B = coordinates_rotate(test_args.C, theta, test_args.C.coord_system.k)
    transformed_field = express(
        test_args.C.coord_system.x * test_args.C.coord_system.i +
        test_args.C.coord_system.y * test_args.C.coord_system.j, B.coord_system)
    with raises(TypeError):
        field_from_sympy_vector(transformed_field, B)


# Test field.apply_to_basis()


def test_basic_apply_to_basis(test_args):
    field = VectorField(test_args.C, lambda p: [p.y, p.x, 0])
    field_space = field.apply_to_basis()
    assert field_space.coordinate_system == test_args.C
    assert field_space.components == [test_args.C.coord_system.y, test_args.C.coord_system.x, 0]


def test_custom_names_apply_to_basis():
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    field = VectorField(C1, lambda p: [p.x, p.y, p.z])
    field_space = field.apply_to_basis()
    assert field_space.coordinate_system == C1
    assert field_space.components == [C1.coord_system.r, C1.coord_system.theta, C1.coord_system.z]


# Test sympy_vector_from_field()


def test_basic_sympy_vector_from_field(test_args):
    field = VectorField(test_args.C, lambda p: [p.y, p.x, 0])
    vector = sympy_vector_from_field(field)
    assert vector == test_args.C.coord_system.y * test_args.C.coord_system.i + test_args.C.coord_system.x * test_args.C.coord_system.j


def test_plain_sympy_vector_from_field(test_args):
    field = VectorField(test_args.C, [1, 2, 0])
    vector = sympy_vector_from_field(field)
    assert vector == 1 * test_args.C.coord_system.i + 2 * test_args.C.coord_system.j


# Test field.apply()


# Result is a function that returns a vector at any point of the trajectory.
# Result is stored in array, where first component of the array is magnitude of the resulting
# vector along X-axis (also called i-vector), second component is magnitude of the resulting
# vector along Y-axis (also called j-vector).
# Input field has X and Y swapped, so as we are moving along X-axis of the trajectory,
# resulting Y component of the vector grows.
def test_basic_field_apply(test_args):
    field = VectorField(test_args.C, lambda p: [p.y, p.x])
    # represents surface
    trajectory = [test_args.C.coord_system.x, test_args.C.coord_system.y]
    trajectory_vectors = field.apply(trajectory)
    assert trajectory_vectors.coordinate_system == test_args.C
    assert trajectory_vectors.components == [test_args.C.coord_system.y, test_args.C.coord_system.x]


def test_parametrized_field_apply(test_args):
    field = field_from_sympy_vector(
        -test_args.C.coord_system.y * test_args.C.coord_system.i +
        test_args.C.coord_system.x * test_args.C.coord_system.j, test_args.C)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_vectors = field.apply(trajectory)
    assert trajectory_vectors.coordinate_system == test_args.C
    assert trajectory_vectors.components == [-1 * parameter, parameter, 0]


def test_sympy_field_apply(test_args):
    field = field_from_sympy_vector(
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
def test_different_coord_systems_field_apply(test_args):
    result_field = VectorField(test_args.C, lambda p: [p.y * p.x, 0, 0])
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    trajectory = [C1.coord_system.r, C1.coord_system.theta]
    trajectory_vectors = result_field.apply(trajectory)
    assert trajectory_vectors.components == [C1.coord_system.r * C1.coord_system.theta, 0, 0]


# Test field_rebase()


def test_basic_field_rebase(test_args):
    field = VectorField(test_args.C, lambda p: [p.x + p.y, 0, 0])
    assert field.coordinate_system == test_args.C
    point = [1, 2, 3]
    point_vector = field.apply(point)
    assert point_vector.components == [3, 0, 0]
    assert point_vector.coordinate_system == test_args.C

    # B is located at [1, 2, 0] origin instead of [0, 0, 0] of test_args.C
    Bi = test_args.C.coord_system.locate_new(
        "B", test_args.C.coord_system.i + 2 * test_args.C.coord_system.j)
    B = CoordinateSystem(test_args.C.coord_system_type, Bi)
    field_rebased = field_rebase(field, B)
    assert field_rebased.basis == [B.coord_system.x, B.coord_system.y, B.coord_system.z]
    assert field_rebased.coordinate_system == B

    transformed_point_vector = field_rebased.apply(point)
    assert transformed_point_vector.components != point_vector.components
    # After rebase field was extended to 3D space
    assert transformed_point_vector.components == [6, 0, 0]


# VectorField invariant does not hold, when applied to some fixed point in space. Use
# 'field_rebase' to let VectorField know about new coordinate system.
def test_invariant_field_rebase_and_apply(test_args):
    field = VectorField(test_args.C, lambda p: [p.x**2 + 2 * p.y**2, 0, 0])
    assert field.coordinate_system == test_args.C
    point = [1, 2]
    p1 = test_args.C.coord_system.origin.locate_new(
        "p1", point[0] * test_args.C.coord_system.i + point[1] * test_args.C.coord_system.j)
    p1_coordinates = p1.express_coordinates(test_args.C.coord_system)
    assert p1_coordinates[0] == point[0]
    assert p1_coordinates[1] == point[1]

    point_vector = field.apply(point)
    assert point_vector.components == [9, 0, 0]

    B = coordinates_rotate(test_args.C, pi / 4, test_args.C.coord_system.k)
    p1_coordinates_in_b = p1.express_coordinates(B.coord_system)
    assert p1_coordinates_in_b[0] != point[0]

    transformed_point = [p1_coordinates_in_b[0], p1_coordinates_in_b[1]]
    transformed_point_vector = field.apply(transformed_point)
    # invariant does not hold if field is not rebased to new coordinate system
    assert transformed_point_vector.components != point_vector.components

    field_rebased = field_rebase(field, B)
    assert field_rebased.coordinate_system == B
    transformed_point_vector = field_rebased.apply(transformed_point)
    # here vector is the same as in 'test_args.C' coordinate system, but it is
    # rotated from the point of view of 'B' coordinate system.
    assert transformed_point_vector.components == [9 * sqrt(2) / 2, -9 * sqrt(2) / 2, 0]


# Test non-cartesian coordinate systems


def test_cylindrical_field_create(test_args):
    field = VectorField(test_args.C, lambda p: [p.x, p.y, 0])
    point = [1, 2]
    point_vector = field.apply(point)
    assert point_vector.components == [1, 2, 0]

    B = coordinates_transform(test_args.C, CoordinateSystem.System.CYLINDRICAL)
    field_rebased = field_rebase(field, B)
    assert field_rebased.coordinate_system == B

    # point should have r = sqrt(5) in polar coordinates
    # theta angle is atan(2/1)
    point_polar = [sqrt(5), atan(2)]
    point_polar_vector = field_rebased.apply(point_polar)
    assert point_polar_vector.components == [sqrt(5), atan(2), 0]

    # now rebase polar vector back to cartesian coordinates and confirm
    # it is the same as original vector
    vector_rebased = vector_rebase(point_polar_vector, test_args.C)
    assert vector_rebased.components == [1, 2, 0]
