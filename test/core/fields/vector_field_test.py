from collections import namedtuple
from typing import Any, List
from pytest import fixture, raises
from sympy import cos, pi, sin, sqrt, symbols
from sympy.vector import CoordSys3D, express, Vector
from test.test_decorators import unsupported_usage
from symplyphysics.core.fields.field_point import FieldPoint
from symplyphysics.core.fields.vector_field import VectorField, field_from_sympy_vector, field_rebase


def _assert_callable(field_: VectorField, size_: int):
    for i in range(size_):
        assert callable(field_.component(i))

def _assert_point(field_: VectorField, point_: FieldPoint, expected_: List[Any]):
    for idx, c in enumerate(field_.components):
        value = c(point_) if callable(c) else c
        assert value == expected_[idx]

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
    _assert_point(field, field_point, [2, 1])
    assert field.basis == []
    assert field.coord_system is None

def test_empty_field():
    field = VectorField()
    assert len(list(field.components)) == 0
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [])

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
    _assert_point(field, field_point, ["string", 1])

@unsupported_usage
def test_wrong_type_value_field():
    field = VectorField("string", lambda p: p.x)
    assert len(list(field.components)) == 2
    assert not callable(field.component(0))
    assert callable(field.component(1))
    field_point = FieldPoint(1, 2, 3)
    # non expression in a field is not processed and returns as is
    _assert_point(field, field_point, ["string", 1])

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
    _assert_point(field, field_point, ["2", 1])

def test_coord_system_field(test_args):
    field = VectorField(lambda p: p.y * p.z, 0, 0, test_args.C)
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [6])
    assert field.basis == [test_args.C.x, test_args.C.y, test_args.C.z]
    assert field.coord_system == test_args.C

# Test field_from_sympy_vector()

def test_basic_vector_to_field_conversion(test_args):
    field = field_from_sympy_vector(test_args.C.x * test_args.C.i + test_args.C.y * test_args.C.j)
    _assert_callable(field, 2)
    assert field.component(2) == 0
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [1, 2, 0])
    assert field.basis == [test_args.C.x, test_args.C.y, test_args.C.z]
    assert field.coord_system == test_args.C

def test_skip_dimension_vector_to_field_conversion(test_args):
    field = field_from_sympy_vector(1 * test_args.C.i + 2 * test_args.C.k)
    assert len(list(field.components)) == 3
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [1, 0, 2])

def test_empty_vector_to_field_conversion():
    field = field_from_sympy_vector(0)
    assert len(list(field.components)) == 0
    # applying empty field to a point results in all zeroes
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [])
    assert field.basis == []

def test_only_integer_vector_to_field_conversion():
    field = field_from_sympy_vector(1)
    # no coordinate system is available from this SymPy vector, so resulting field is empty
    assert len(list(field.components)) == 0
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [])
    assert field.basis == []

def test_only_scalar_to_field_conversion(test_args):
    with raises(TypeError):
        field_from_sympy_vector(test_args.C.x)

# different coordinate systems in parameters are not supported
def test_different_coord_systems_vector_to_field_conversion(test_args):
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    with raises(TypeError):
        field_from_sympy_vector(test_args.C.x * test_args.C.i + 2 * C1.phi * C1.j)
    with raises(TypeError):
        field_from_sympy_vector(test_args.C.x * C1.i)

def test_custom_names_vector_to_field_conversion():
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    field = field_from_sympy_vector(C1.r * C1.i + 2 * C1.phi * C1.j)
    _assert_callable(field, 2)
    assert field.component(2) == 0
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [1, 4, 0])
    assert field.basis == [C1.r, C1.phi, C1.z]
    assert field.coord_system == C1

def test_rotate_coordinates_vector_to_field_conversion(test_args):
    sympy_vector_field = test_args.C.x * test_args.C.i + test_args.C.y * test_args.C.j
    field = field_from_sympy_vector(sympy_vector_field)
    _assert_callable(field, 2)
    assert field.component(2) == 0
    field_point = FieldPoint(1, 2, 3)
    _assert_point(field, field_point, [1, 2, 0])
    theta = symbols("theta")
    B = test_args.C.orient_new_axis('B', theta, test_args.C.k)
    transformed_vector = express(sympy_vector_field, B, variables=True)
    result_transformed_field = field_from_sympy_vector(transformed_vector)
    _assert_point(result_transformed_field, field_point,
        [(sin(theta) + 2 * cos(theta)) * sin(theta) + (cos(theta) - 2 * sin(theta)) * cos(theta),
        (sin(theta) + 2 * cos(theta)) * cos(theta) - (cos(theta) - 2 * sin(theta)) * sin(theta),
        0])
    assert result_transformed_field.basis == [B.x, B.y, B.z]
    assert result_transformed_field.coord_system == B

# when we express SymPy vector field to another coordinate system and base scalars (C.x, C.y) are
# not transformed, we get multiple coordinate systems, which is unsupported by sympy_vector_to_field
def test_rotate_coordinates_without_variables_vector_to_field_conversion(test_args):
    theta = symbols("theta")
    B = test_args.C.orient_new_axis('B', theta, test_args.C.k)
    transformed_field = express(test_args.C.x * test_args.C.i + test_args.C.y * test_args.C.j, B)
    with raises(TypeError):
        field_from_sympy_vector(transformed_field)

# Test field.apply_to_basis()

def test_basic_apply_to_basis(test_args):
    field = VectorField(lambda p: p.y, lambda p: p.x, 0, test_args.C)
    field_space = field.apply_to_basis()
    assert field_space == [test_args.C.y, test_args.C.x]

def test_custom_names_apply_to_basis():
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    field = VectorField(lambda p: p.x, lambda p: p.y, lambda p: p.z, C1)
    field_space = field.apply_to_basis()
    assert field_space == [C1.r, C1.phi, C1.z]

def test_spherical_apply_to_basis():
    C1 = CoordSys3D("C1", transformation="spherical")
    field = VectorField(lambda p: p.x, lambda p: p.y, lambda p: p.z, C1)
    field_space = field.apply_to_basis()
    assert field_space == [C1.r, C1.theta, C1.phi]

def test_empty_basis_apply_to_basis():
    field = VectorField(lambda p: p.y * p.x)
    field_space = field.apply_to_basis()
    assert field_space == [0]

# Test field.apply()

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
    trajectory_vectors = field.apply(trajectory)
    assert len(trajectory_vectors) == 2
    assert trajectory_vectors == [test_args.C.y, test_args.C.x]

# Coordinate system is not necessary to apply field.
def test_parametrized_no_coord_system_field_apply():
    field = VectorField(lambda p: -p.y, lambda p: p.x)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_vectors = field.apply(trajectory)
    assert len(trajectory_vectors) == 2
    assert trajectory_vectors == [-parameter, parameter]

def test_parametrized_field_apply(test_args):
    field = field_from_sympy_vector(-test_args.C.y * test_args.C.i + test_args.C.x * test_args.C.j)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_vectors = field.apply(trajectory)
    assert len(trajectory_vectors) == 3
    assert trajectory_vectors == [-parameter, parameter, 0]

def test_sympy_field_apply(test_args):
    field = field_from_sympy_vector(-test_args.C.y * test_args.C.i + test_args.C.x * test_args.C.j)
    trajectory = [test_args.C.x, test_args.C.y]
    trajectory_vectors = field.apply(trajectory)
    assert len(trajectory_vectors) == 3
    assert trajectory_vectors == [-test_args.C.y, test_args.C.x, 0]

# VectorField checks that coordinate systems of field and trajectory should be same, if they
# are set. VectorField is not rebased automatically and should be rebased to the same coordinate
# system as in trajectory with 'field_rebase'.
def test_different_coord_systems_field_apply(test_args):
    result_field = VectorField(lambda p: p.y * p.x, 0, 0, test_args.C)
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    trajectory = [C1.r, C1.phi]
    with raises(TypeError):
        result_field.apply(trajectory)

# Test field_rebase()

def test_basic_field_rebase(test_args):
    field = VectorField(lambda p: p.x + p.y, 0, 0, test_args.C)
    point = [1, 2, 3]
    point_vector = field.apply(point)
    assert point_vector == [3]
    assert field.coord_system == test_args.C

    # B is located at [1, 2, 0] origin instead of [0, 0, 0] of test_args.C
    B = test_args.C.locate_new('B', test_args.C.i + 2 * test_args.C.j)
    field_rebased = field_rebase(field, B)
    assert field_rebased.basis == [B.x, B.y, B.z]
    assert field_rebased.coord_system == B

    transformed_point_vector = field_rebased.apply(point)
    assert transformed_point_vector != point_vector
    # After rebase field was extended to 3D space
    assert transformed_point_vector == [6, 0, 0]

# VectorField invariant does not hold, when applied to some fixed point in space. Use
# 'field_rebase' to let VectorField know about new coordinate system.
def test_invariant_field_rebase_and_apply(test_args):
    field = VectorField(lambda p: p.x**2 + 2 * p.y**2, 0, 0, test_args.C)
    point = [1, 2]
    p1 = test_args.C.origin.locate_new('p1', point[0] * test_args.C.i + point[1] * test_args.C.j)
    p1_coordinates = p1.express_coordinates(test_args.C)
    assert p1_coordinates[0] == point[0]
    assert p1_coordinates[1] == point[1]

    point_vector = field.apply(point)
    assert point_vector == [9]
    assert field.coord_system == test_args.C

    B = test_args.C.orient_new_axis('B', pi/4, test_args.C.k)
    p1_coordinates_in_b = p1.express_coordinates(B)
    assert p1_coordinates_in_b[0] != point[0]

    transformed_point = [ p1_coordinates_in_b[0], p1_coordinates_in_b[1] ]
    transformed_point_vector = field.apply(transformed_point)
    # invariant does not hold if field is not rebased to new coordinate system
    assert transformed_point_vector != point_vector

    field_rebased = field_rebase(field, B)
    assert field_rebased.coord_system == B
    transformed_point_vector = field_rebased.apply(transformed_point)
    # here vector is the same as in 'test_args.C' coordinate system, but it is
    # rotated from the point of view of 'B' coordinate system.
    assert transformed_point_vector == [9 * sqrt(2) / 2, -9 * sqrt(2) / 2, 0]

# Field is not rebased if no original coordinate system was set.
def test_no_coord_system_field_rebase(test_args):
    field = VectorField(lambda p: p.x + p.y)
    assert field.coord_system is None
    point = [1, 2]
    point_vector = field.apply(point)
    assert point_vector == [3]

    B = test_args.C.locate_new('B', test_args.C.i + 2 * test_args.C.j)
    field_rebased = field_rebase(field, B)
    assert field_rebased.basis == [B.x, B.y, B.z]
    assert field_rebased.coord_system == B

    point_vector = field_rebased.apply(point)
    assert point_vector == [3]

# Field is not rebased if no target coordinate system was set.
def test_no_target_coord_system_field_rebase(test_args):
    field = VectorField(lambda p: p.x + p.y, 0, 0, test_args.C)
    assert field.coord_system == test_args.C
    point = [1, 2]
    point_vector = field.apply(point)
    assert point_vector == [3]

    field_rebased = field_rebase(field, None)
    assert field_rebased.coord_system == None

    point_vector = field_rebased.apply(point)
    assert point_vector == [3]
