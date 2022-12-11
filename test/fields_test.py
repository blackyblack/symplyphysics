from collections import namedtuple
import functools
from typing import Any, List
from pytest import fixture, raises
from sympy.vector import CoordSys3D
from symplyphysics.fields import FieldPoint, VectorField, apply_field, field_from_unit_vector


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

def test_basic_array_to_field_conversion(test_args):
    result_field = field_from_unit_vector([test_args.C.x, test_args.C.y])
    assert len(list(result_field.components)) == 2
    _assert_callable(result_field, 2)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [1, 2, 0])

def test_skip_dimension_array_to_field_conversion():
    result_field = field_from_unit_vector([1, 0, 2])
    assert len(list(result_field.components)) == 3
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [1, 0, 2])

def test_4d_array_to_field_conversion(test_args):
    result_field = field_from_unit_vector([test_args.C.x, 1, 2, 5])
    assert len(list(result_field.components)) == 4
    _assert_callable(result_field, 1)
    field_point = FieldPoint(1, 2, 3)
    field_point.set_coordinate(3, 4)
    _assert_point(result_field, field_point, [1, 1, 2, 5])

def test_empty_array_to_field_conversion():
    result_field = field_from_unit_vector([])
    assert len(list(result_field.components)) == 0
    # applying empty field to a point results in all zeroes
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [0, 0, 0])

# TODO: should fail when different coordinate systems in a field
def test_different_coord_systems_array_to_field_conversion(test_args):
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    result_field = field_from_unit_vector([test_args.C.x, 2 * C1.phi])
    assert len(list(result_field.components)) == 2
    _assert_callable(result_field, 2)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [1, 4, 0])

def test_custom_names_array_to_field_conversion():
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    result_field = field_from_unit_vector([C1.r, 2 * C1.phi])
    assert len(list(result_field.components)) == 2
    _assert_callable(result_field, 2)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [1, 4, 0])

def test_basic_field_apply(test_args):
    field = VectorField(lambda p: p.y, lambda p: p.x)
    # represents surface
    trajectory = [test_args.C.x, test_args.C.y]
    trajectory_vector = apply_field(test_args.C, field, trajectory)
    assert len(trajectory_vector) == 2
    # result vector has X component proportional to Y coordinate of a point,
    # Y component proportional to X coordinate of a point
    assert trajectory_vector == [test_args.C.y, test_args.C.x]
