from collections import namedtuple
from pytest import fixture, raises
from sympy import cos, pi, sin, symbols
from sympy.vector import CoordSys3D, express
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from test.test_decorators import unsupported_usage
from symplyphysics.core.fields.field_point import FieldPoint
from symplyphysics.core.fields.scalar_field import ScalarField, field_from_sympy_vector, field_rebase, sympy_expression_to_field_function


@fixture
def test_args():
    C = CoordinateSystem()
    Args = namedtuple("Args", ["C"])
    return Args(C=C)

# Test sympy_expression_to_field_function()

def test_basic_sympy_expression_to_field_function(test_args):
    field_function = sympy_expression_to_field_function(test_args.C.coord_system.x + test_args.C.coord_system.y, test_args.C)
    assert callable(field_function)
    field_point = FieldPoint(1, 2, 3)
    assert field_function(field_point) == 3

def test_empty_sympy_expression_to_field_function():
    field_function = sympy_expression_to_field_function(0)
    assert field_function == 0

# Any non SymPy expression is returned without modification
def test_integer_sympy_expression_to_field_function():
    field_function = sympy_expression_to_field_function(1)
    assert field_function == 1

# Only scalars from requested coordinate system are being applied
def test_partially_different_coord_systems_sympy_expression_to_field_function(test_args):
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    field_function = sympy_expression_to_field_function(test_args.C.coord_system.x + 2 * C1.phi, test_args.C)
    assert callable(field_function)
    field_point = FieldPoint(1, 2, 3)
    assert field_function(field_point) == 1 + 2 * C1.phi

# Test ScalarField constructor

def test_basic_field():
    field = ScalarField(lambda p: p.y * p.z)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 6
    assert field.basis == []
    assert field.coordinate_system is None
    assert len(field.components) == 1
    assert callable(field.components[0])

def test_empty_field():
    field = ScalarField()
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 0
    # Scalar field components size is always 1
    assert len(field.components) == 1
    assert field.components[0] == 0

def test_4d_point_field():
    field = ScalarField(lambda p: p.coordinate(3))
    field_point = FieldPoint(1, 2, 3)
    field_point.set_coordinate(3, 4)
    assert field(field_point) == 4

@unsupported_usage
def test_wrong_type_lambda_field():
    field = ScalarField(lambda p: "string")
    field_point = FieldPoint(1, 2, 3)
    # non expression lambda in a field is not processed and returns as is
    assert field(field_point) == "string"

@unsupported_usage
def test_wrong_type_value_field():
    field = ScalarField("string")
    field_point = FieldPoint(1, 2, 3)
    # non expression in a field is not processed and returns as is
    assert field(field_point) == "string"

@unsupported_usage
def test_invalid_lambda_field():
    field = ScalarField(lambda p: p.y + "string")
    field_point = FieldPoint(1, 2, 3)
    # cannot add integer and string in field lambda
    with raises(TypeError):
        field(field_point)

@unsupported_usage
def test_effect_in_lambda_field():
    field = ScalarField(lambda p: "{}".format(p.y))
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == "2"

def test_coord_system_field(test_args):
    field = ScalarField(lambda p: p.y * p.z, test_args.C)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 6
    assert field.basis == [test_args.C.coord_system.x, test_args.C.coord_system.y, test_args.C.coord_system.z]
    assert field.coordinate_system == test_args.C

# Test field_from_sympy_vector()

def test_basic_vector_to_field_conversion(test_args):
    field = field_from_sympy_vector(test_args.C.coord_system.x + test_args.C.coord_system.y, test_args.C)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 3
    assert field.basis == [test_args.C.coord_system.x, test_args.C.coord_system.y, test_args.C.coord_system.z]
    assert field.coordinate_system == test_args.C

# Coordinate system base vectors (C.i, C.j) are not being processed by ScalarField,
# so as a result of applying field they are kept untouched, like all other free variables in SymPy expression.
def test_dimensional_vector_to_field_conversion(test_args):
    field = field_from_sympy_vector(test_args.C.coord_system.x * test_args.C.coord_system.i + test_args.C.coord_system.y * test_args.C.coord_system.j, test_args.C)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == test_args.C.coord_system.i + 2 * test_args.C.coord_system.j

def test_empty_vector_to_field_conversion():
    field = field_from_sympy_vector(0)
    # applying empty field to a point results in zero value
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 0
    assert field.basis == []

def test_only_integer_vector_to_field_conversion():
    field = field_from_sympy_vector(1)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 1
    assert field.basis == []

def test_custom_names_vector_to_field_conversion():
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    field = field_from_sympy_vector(C1.coord_system.r + 2 * C1.coord_system.theta, C1)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 5
    assert field.basis == [C1.coord_system.r, C1.coord_system.theta, C1.coord_system.z]
    assert field.coordinate_system == C1

def test_rotate_coordinates_vector_to_field_conversion(test_args):
    sympy_vector_field = test_args.C.coord_system.x + test_args.C.coord_system.y
    field = field_from_sympy_vector(sympy_vector_field, test_args.C)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 3
    theta = symbols("theta")
    B_inner = test_args.C.coord_system.orient_new_axis('B', theta, test_args.C.coord_system.k)
    B = CoordinateSystem(test_args.C.coord_system_type, B_inner)
    transformed_vector = express(sympy_vector_field, B_inner, variables=True)
    result_transformed_field = field_from_sympy_vector(transformed_vector, B)
    assert transformed_vector == B.coord_system.x * sin(theta) + B.coord_system.x * cos(theta) - B.coord_system.y * sin(theta) + B.coord_system.y * cos(theta)
    assert result_transformed_field(field_point) == sin(theta) + cos(theta) - 2 * sin(theta) + 2 * cos(theta)
    assert result_transformed_field.basis == [B.coord_system.x, B.coord_system.y, B.coord_system.z]
    assert result_transformed_field.coordinate_system == B

# when we express SymPy Vector to another coordinate system and base scalars (C.x, C.y) are
# not transformed, we have same base scalars as before.
def test_rotate_coordinates_without_variables_vector_to_field_conversion(test_args):
    theta = symbols("theta")
    B_inner = test_args.C.coord_system.orient_new_axis('B', theta, test_args.C.coord_system.k)
    B = CoordinateSystem(test_args.C.coord_system_type, B_inner)
    transformed_field = express(test_args.C.coord_system.x + test_args.C.coord_system.y, B_inner)
    field = field_from_sympy_vector(transformed_field, B)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == test_args.C.coord_system.x + test_args.C.coord_system.y

# Test field.apply()

# Result is a function that returns a scalar value at any point of the trajectory.
# Input field has X * Y, so as we are moving along X-axis or Y-axis of the trajectory,
# resulting value grows proportionally to X * Y.
def test_basic_field_apply(test_args):
    field = ScalarField(lambda p: p.y * p.x)
    # represents surface
    trajectory = [test_args.C.coord_system.x, test_args.C.coord_system.y]
    trajectory_value = field.apply(trajectory)
    assert trajectory_value == test_args.C.coord_system.y * test_args.C.coord_system.x

# Coordinate system is not necessary to apply field.
def test_parametrized_no_coord_system_field_apply():
    field = ScalarField(lambda p: -p.y)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_value = field.apply(trajectory)
    assert trajectory_value == -parameter

def test_parametrized_field_apply(test_args):
    result_field = field_from_sympy_vector(-test_args.C.coord_system.y + 2 * test_args.C.coord_system.x, test_args.C)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_value = result_field.apply(trajectory)
    assert trajectory_value == parameter

def test_sympy_field_apply(test_args):
    result_field = field_from_sympy_vector(-test_args.C.coord_system.y + test_args.C.coord_system.x, test_args.C)
    trajectory = [test_args.C.coord_system.x, test_args.C.coord_system.y]
    trajectory_value = result_field.apply(trajectory)
    assert trajectory_value == -test_args.C.coord_system.y + test_args.C.coord_system.x

def test_uncallable_field_apply(test_args):
    result_field = ScalarField(1)
    trajectory = [test_args.C.coord_system.x, test_args.C.coord_system.y]
    trajectory_value = result_field.apply(trajectory)
    assert trajectory_value == 1

# ScalarField checks that coordinate systems of field and trajectory should be same, if they
# are set. ScalarField is not rebased automatically and should be rebased to the same coordinate
# system as in trajectory with 'field_rebase'.
def test_different_coord_systems_field_apply(test_args):
    result_field = ScalarField(lambda p: p.y * p.x, test_args.C)
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    trajectory = [C1.r, C1.phi]
    with raises(TypeError):
        result_field.apply(trajectory)

# While ScalarField should contain information about coordinate system and
# can be rebased to new coordinate system with ScalarField invariance, it is not
# necessary to do so. Instead user can define a trajectory, that contains this point, transform
# trajectory to new coordinate system and apply field to it.
# ScalarField invariant will hold.
def test_invariant_transformed_trajectory_field_apply(test_args):
    field = ScalarField(lambda p: p.x**2 + 2 * p.y**2)
    point = [1, 2]
    trajectory = [test_args.C.x, test_args.C.y + 5]
    trajectory_value = field.apply(trajectory)
    assert trajectory_value == test_args.C.x**2 + 2 * (test_args.C.y + 5)**2
    assert trajectory_value.subs({test_args.C.x: point[0], test_args.C.y: point[1]}) == 99

    B = test_args.C.orient_new_axis('B', pi/4, test_args.C.k)
    transformed_trajectory = [
        express(trajectory[0], B, variables=True),
        express(trajectory[1], B, variables=True)]
    transformed_trajectory_value = field.apply(transformed_trajectory)

    p1 = test_args.C.origin.locate_new('p1', point[0] * test_args.C.i + point[1] * test_args.C.j)
    p1_coordinates = p1.express_coordinates(test_args.C)
    assert p1_coordinates[0] == point[0]
    assert p1_coordinates[1] == point[1]

    p1_coordinates_in_b = p1.express_coordinates(B)
    assert p1_coordinates_in_b[0] != point[0]

    assert transformed_trajectory_value.subs({B.x: p1_coordinates_in_b[0], B.y: p1_coordinates_in_b[1]}) == 99

# Test field.apply_to_basis()

def test_basic_field_apply_to_basis(test_args):
    field = ScalarField(lambda p: p.y * p.x, test_args.C)
    field_space = field.apply_to_basis()
    assert field_space == test_args.C.y * test_args.C.x

def test_empty_basis_apply_to_basis():
    field = ScalarField(lambda p: p.y * p.x)
    field_space = field.apply_to_basis()
    assert field_space == 0

# Test field_rebase()

def test_basic_field_rebase(test_args):
    field = ScalarField(lambda p: p.x + p.y, test_args.C)
    point = [1, 2]
    point_value = field.apply(point)
    assert point_value == 3
    assert field.coord_system == test_args.C

    # B is located at [1, 2] origin instead of [0, 0] of test_args.C
    B = test_args.C.locate_new('B', test_args.C.i + 2 * test_args.C.j)
    field_rebased = field_rebase(field, B)
    assert field_rebased.basis == [B.x, B.y, B.z]
    assert field_rebased.coord_system == B
    # Original field is not changed
    assert field.basis == [test_args.C.x, test_args.C.y, test_args.C.z]
    assert field.coord_system == test_args.C

    transformed_point_value = field_rebased.apply(point)
    assert transformed_point_value != point_value
    assert transformed_point_value == 6

# ScalarField invariant does not hold, when applied to some fixed point in space. Use
# 'field_rebase' to let ScalarField know about new coordinate system.
def test_invariant_field_rebase_and_apply(test_args):
    field = ScalarField(lambda p: p.x**2 + 2 * p.y**2, test_args.C)
    point = [1, 2]
    p1 = test_args.C.origin.locate_new('p1', point[0] * test_args.C.i + point[1] * test_args.C.j)
    p1_coordinates = p1.express_coordinates(test_args.C)
    assert p1_coordinates[0] == point[0]
    assert p1_coordinates[1] == point[1]

    point_value = field.apply(point)
    assert point_value == 9

    B = test_args.C.orient_new_axis('B', pi/4, test_args.C.k)
    p1_coordinates_in_b = p1.express_coordinates(B)
    assert p1_coordinates_in_b[0] != point[0]

    transformed_point = [ p1_coordinates_in_b[0], p1_coordinates_in_b[1] ]
    transformed_point_value = field.apply(transformed_point)
    # invariant does not hold if field is not rebased to new coordinate system
    assert transformed_point_value != point_value

    field_rebased = field_rebase(field, B)
    transformed_point_value = field_rebased.apply(transformed_point)
    assert transformed_point_value == point_value

# Field is not rebased if no original coordinate system was set.
def test_no_coord_system_field_rebase(test_args):
    field = ScalarField(lambda p: p.x + p.y)
    assert field.coord_system is None
    point = [1, 2]
    point_value = field.apply(point)
    assert point_value == 3

    B = test_args.C.locate_new('B', test_args.C.i + 2 * test_args.C.j)
    field_rebased = field_rebase(field, B)
    assert field_rebased.basis == [B.x, B.y, B.z]
    assert field_rebased.coord_system == B

    point_value = field_rebased.apply(point)
    assert point_value == 3

# Field is not rebased if no target coordinate system was set.
def test_no_target_coord_system_field_rebase(test_args):
    field = ScalarField(lambda p: p.x + p.y, test_args.C)
    assert field.coord_system == test_args.C
    point = [1, 2]
    point_value = field.apply(point)
    assert point_value == 3

    field_rebased = field_rebase(field, None)
    assert field_rebased.coord_system == None

    point_value = field_rebased.apply(point)
    assert point_value == 3
