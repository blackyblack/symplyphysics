from collections import namedtuple
from typing import Any, List
from pytest import fixture
from sympy.vector import CoordSys3D
from symplyphysics.fields import FieldPoint, VectorField, field_from_unit_vector


@fixture
def test_args():
    C = CoordSys3D("C")
    Args = namedtuple("Args", ["C"])
    return Args(C=C)


def _assert_callable(field_: VectorField, size_: int):
    for i in range(size_):
        assert callable(field_.component(i))

def _assert_point(field_: VectorField, point_: FieldPoint, expected_: List[Any]):
    for idx, c in enumerate(field_.components):
        value = c(point_) if callable(c) else c
        assert value == expected_[idx]

def test_basic_array_to_field_conversion(test_args):
    result_field = field_from_unit_vector(test_args.C, [test_args.C.x, test_args.C.y])
    _assert_callable(result_field, 2)
    assert result_field.component(2) == 0
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [1, 2, 0])

def test_skip_dimension_array_to_field_conversion(test_args):
    result_field = field_from_unit_vector(test_args.C, [1, 0, 2])
    _assert_callable(result_field, 3)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(result_field, field_point, [1, 0, 2])

def test_4d_array_to_field_conversion(test_args):
    result_field = field_from_unit_vector(test_args.C, [test_args.C.x, 1, 2, 5])
    # this coordinate system only has 3 dimensions
    assert len(list(result_field.components)) == 3
    _assert_callable(result_field, 3)
    field_point = FieldPoint(1, 2, 3)
    # add 4th coordinate and see it does not affect anything
    field_point.set_coordinate(3, 4)
    _assert_point(result_field, field_point, [1, 1, 2])
