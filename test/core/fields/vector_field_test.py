from collections import namedtuple
import functools
from typing import Any, List
from pytest import fixture, raises
from sympy import cos, sin, symbols
from sympy.vector import CoordSys3D, express
from symplyphysics.core.fields.field_point import FieldPoint
from symplyphysics.core.fields.vector_field import VectorField, apply_field, apply_field_to_coord_system, sympy_vector_to_field


def _assert_callable(field_: VectorField, size_: int):
    for i in range(size_):
        assert callable(field_.component(i))

def _assert_point(field_: VectorField, point_: FieldPoint, expected_: List[Any]):
    for idx, c in enumerate(field_.components):
        value = c(point_) if callable(c) else c
        assert value == expected_[idx]

# This decorator is used to mark tests that show unsupported usage of VectorField.
# You're not supposed to use VectorField this way, but we cannot enforce it by Python code.
# Use tests marked with `unsupported_usage` as a reference on how field behaves when passed incorrect input.
def unsupported_usage(func):
    @functools.wraps(func)
    def wrapper_(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper_


@fixture
def test_args():
    C = CoordSys3D("C")
    Args = namedtuple("Args", ["C"])
    return Args(C=C)

# Test VectorField constructor

def test_basic_field():
    field = VectorField(lambda p: p.y, lambda p: p.x)
    assert len(list(field.components)) == 2
    _assert_callable(field, 2)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [2, 1, 0])

def test_empty_field():
    field = VectorField()
    assert len(list(field.components)) == 0
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [0, 0, 0])

def test_4d_field():
    field = VectorField(lambda p: p.x, lambda p: p.y, lambda p: p.z)
    field.set_component(3, lambda p: p.x)
    assert len(list(field.components)) == 4
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [1, 2, 3, 1])

def test_4d_point_field():
    field = VectorField(lambda p: p.x, lambda p: p.y, lambda p: p.z)
    field.set_component(3, lambda p: p.coordinate(3))
    assert len(list(field.components)) == 4
    field_point = FieldPoint(1, 2, 3)
    field_point.set_coordinate(3, 4)
    _assert_point(field, field_point, [1, 2, 3, 4])

@unsupported_usage
def test_wrong_type_lambda_field():
    field = VectorField(lambda p: "string", lambda p: p.x)
    assert len(list(field.components)) == 2
    _assert_callable(field, 2)
    field_point = FieldPoint(1, 2, 3)
    # non expression lambda in a field is not processed and returns as is
    _assert_point(field, field_point, ["string", 1, 0])

@unsupported_usage
def test_wrong_type_value_field():
    field = VectorField("string", lambda p: p.x)
    assert len(list(field.components)) == 2
    assert not callable(field.component(0))
    assert callable(field.component(1))
    field_point = FieldPoint(1, 2, 3)
    # non expression in a field is not processed and returns as is
    _assert_point(field, field_point, ["string", 1, 0])

@unsupported_usage
def test_invalid_lambda_field():
    field = VectorField(lambda p: p.y + "string", lambda p: p.x)
    assert len(list(field.components)) == 2
    _assert_callable(field, 2)
    field_point = FieldPoint(1, 2, 3)
    # cannot add integer and string in field lambda
    with raises(TypeError):
        field.component(0)(field_point)

@unsupported_usage
def test_effect_in_lambda_field():
    field = VectorField(lambda p: "{}".format(p.y), lambda p: p.x)
    assert len(list(field.components)) == 2
    _assert_callable(field, 2)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, ["2", 1, 0])

# Test sympy_vector_to_field()

def test_basic_vector_to_field_conversion(test_args):
    result_field = sympy_vector_to_field(test_args.C.x * test_args.C.i + test_args.C.y * test_args.C.j)
    _assert_callable(result_field, 3)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [1, 2, 0])

def test_skip_dimension_vector_to_field_conversion(test_args):
    result_field = sympy_vector_to_field(1 * test_args.C.i + 2 * test_args.C.k)
    assert len(list(result_field.components)) == 3
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [1, 0, 2])

def test_empty_vector_to_field_conversion():
    result_field = sympy_vector_to_field(0)
    assert len(list(result_field.components)) == 0
    # applying empty field to a point results in all zeroes
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [0, 0, 0])

def test_only_integer_vector_to_field_conversion():
    result_field = sympy_vector_to_field(1)
    # no coordinate system is available from this SymPy vector, so resulting field is empty
    assert len(list(result_field.components)) == 0
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [0, 0, 0])

def test_only_scalar_to_field_conversion(test_args):
    with raises(TypeError):
        sympy_vector_to_field(test_args.C.x)

# different coordinate systems in parameters are not supported
def test_different_coord_systems_vector_to_field_conversion(test_args):
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    with raises(TypeError):
        sympy_vector_to_field(test_args.C.x * test_args.C.i + 2 * C1.phi * C1.j)
    with raises(TypeError):
        sympy_vector_to_field(test_args.C.x * C1.i)

def test_custom_names_vector_to_field_conversion():
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    result_field = sympy_vector_to_field(C1.r * C1.i + 2 * C1.phi * C1.j)
    _assert_callable(result_field, 3)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [1, 4, 0])

def test_rotate_coordinates_vector_to_field_conversion(test_args):
    sympy_vector_field = test_args.C.x * test_args.C.i + test_args.C.y * test_args.C.j
    result_field = sympy_vector_to_field(sympy_vector_field)
    _assert_callable(result_field, 3)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [1, 2, 0])
    theta = symbols("theta")
    B = test_args.C.orient_new_axis('B', theta, test_args.C.k)
    transformed_vector = express(sympy_vector_field, B, variables=True)
    result_transformed_field = sympy_vector_to_field(transformed_vector)
    _assert_callable(result_transformed_field, 3)

    _assert_point(result_transformed_field, field_point,
        [(sin(theta) + 2 * cos(theta)) * sin(theta) + (cos(theta) - 2 * sin(theta)) * cos(theta),
        (sin(theta) + 2 * cos(theta)) * cos(theta) - (cos(theta) - 2 * sin(theta)) * sin(theta),
        0])

# when we express sympy vector field to another coordinate system and base scalars (C.x, C.y) are
# not transformed, we get multiple coordinate systems, which is unsupported by sympy_vector_to_field
def test_rotate_coordinates_without_variables_vector_to_field_conversion(test_args):
    theta = symbols("theta")
    B = test_args.C.orient_new_axis('B', theta, test_args.C.k)
    transformed_field = express(test_args.C.x * test_args.C.i + test_args.C.y * test_args.C.j, B)
    with raises(TypeError):
        sympy_vector_to_field(transformed_field)

# Test apply_field_to_coord_system()

def test_basic_apply_field_to_coord_system(test_args):
    field = VectorField(lambda p: p.y, lambda p: p.x)
    field_space = apply_field_to_coord_system(field, test_args.C)
    assert field_space == [test_args.C.y, test_args.C.x]

def test_custom_names_apply_field_to_coord_system():
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    field = VectorField(lambda p: p.x, lambda p: p.y, lambda p: p.z)
    field_space = apply_field_to_coord_system(field, C1)
    assert field_space == [C1.r, C1.phi, C1.z]

def test_spherical_apply_field_to_coord_system():
    C1 = CoordSys3D("C1", transformation="spherical")
    field = VectorField(lambda p: p.x, lambda p: p.y, lambda p: p.z)
    field_space = apply_field_to_coord_system(field, C1)
    assert field_space == [C1.r, C1.theta, C1.phi]

# Test apply_field()

# Result is a function that returns a vector at any point of the trajectory.
# Result is stored in array, where first component of the array is magnitude of the resulting
# vector along X-axis (also called i-vector), second component is magnitude of the resulting
# vector along Y-axis (also called j-vector).
# Input field has X and Y swapped, so as we are moving along X-axis of the trajectory,
# resulting Y component of the vector grows.
def test_basic_field_apply(test_args):
    field = VectorField(lambda p: p.y, lambda p: p.x)
    # represents surface
    trajectory = [test_args.C.x, test_args.C.y]
    trajectory_vectors = apply_field(field, trajectory)
    assert len(trajectory_vectors) == 2
    assert trajectory_vectors == [test_args.C.y, test_args.C.x]

# Coordinate system is not necessary to apply field.
def test_parametrized_no_coord_system_field_apply():
    field = VectorField(lambda p: -p.y, lambda p: p.x)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_vectors = apply_field(field, trajectory)
    assert len(trajectory_vectors) == 2
    assert trajectory_vectors == [-parameter, parameter]

def test_parametrized_field_apply(test_args):
    result_field = sympy_vector_to_field(-test_args.C.y * test_args.C.i + test_args.C.x * test_args.C.j)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_vectors = apply_field(result_field, trajectory)
    assert len(trajectory_vectors) == 3
    assert trajectory_vectors == [-parameter, parameter, 0]

def test_sympy_field_apply(test_args):
    result_field = sympy_vector_to_field(-test_args.C.y * test_args.C.i + test_args.C.x * test_args.C.j)
    trajectory = [test_args.C.x, test_args.C.y]
    trajectory_vectors = apply_field(result_field, trajectory)
    assert len(trajectory_vectors) == 3
    assert trajectory_vectors == [-test_args.C.y, test_args.C.x, 0]
