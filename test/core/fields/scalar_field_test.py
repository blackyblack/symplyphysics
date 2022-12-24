from collections import namedtuple
from pytest import fixture, raises
from sympy import cos, sin, symbols
from sympy.vector import CoordSys3D, express
from test.test_decorators import unsupported_usage
from symplyphysics.core.fields.field_point import FieldPoint
from symplyphysics.core.fields.scalar_field import ScalarField, apply_field, sympy_expression_to_field_function, sympy_vector_to_field, extract_coord_system_from_sympy_vector


@fixture
def test_args():
    C = CoordSys3D("C")
    Args = namedtuple("Args", ["C"])
    return Args(C=C)

# Test ScalarField constructor

def test_basic_field():
    field = ScalarField(lambda p: p.y * p.z)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 6

def test_empty_field():
    field = ScalarField()
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 0

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

# Test sympy_expression_to_field_function()

def test_basic_sympy_expression_to_field_function(test_args):
    field_function = sympy_expression_to_field_function(test_args.C, test_args.C.x + test_args.C.y)
    assert callable(field_function)
    field_point = FieldPoint(1, 2, 3)
    assert field_function(field_point) == 3

def test_empty_sympy_expression_to_field_function(test_args):
    field_function = sympy_expression_to_field_function(test_args.C, 0)
    assert field_function == 0

# Empty coordinate system results to unchanged original expression
def test_empty_coord_system_sympy_expression_to_field_function(test_args):
    field_function = sympy_expression_to_field_function(None, test_args.C.x + test_args.C.y)
    assert field_function == test_args.C.x + test_args.C.y

# Any non SymPy expression is returned without modification
def test_integer_sympy_expression_to_field_function(test_args):
    field_function = sympy_expression_to_field_function(test_args.C, 1)
    assert field_function == 1

# This should not happen as coordinate system is detected automatically with extract_coord_system_from_sympy_vector(),
# but as public method it can be called with wrong parameters. In this case, only correct coordinates of SymPy expression
# will be modified.
def test_different_coord_systems_sympy_expression_to_field_function(test_args):
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    field_function = sympy_expression_to_field_function(test_args.C, C1.r + 2 * C1.phi)
    field_point = FieldPoint(1, 2, 3)
    assert field_function(field_point) == C1.r + 2 * C1.phi

def test_partially_different_coord_systems_sympy_expression_to_field_function(test_args):
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    field_function = sympy_expression_to_field_function(test_args.C, test_args.C.x + 2 * C1.phi)
    field_point = FieldPoint(1, 2, 3)
    assert field_function(field_point) == 1 + 2 * C1.phi

# Test extract_coord_system_from_sympy_vector()

def test_basic_extract_coord_system_from_sympy_vector(test_args):
    coord_system = extract_coord_system_from_sympy_vector(test_args.C.x + test_args.C.y)
    assert coord_system == test_args.C

def test_empty_extract_coord_system_from_sympy_vector():
    coord_system = extract_coord_system_from_sympy_vector(1)
    assert coord_system is None

def test_different_coord_systems_extract_coord_system_from_sympy_vector(test_args):
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    with raises(TypeError):
        extract_coord_system_from_sympy_vector(test_args.C.x + 2 * C1.phi)
    with raises(TypeError):
        extract_coord_system_from_sympy_vector(test_args.C.x * C1.i)

# Test sympy_vector_to_field()

def test_basic_vector_to_field_conversion(test_args):
    field = sympy_vector_to_field(test_args.C.x + test_args.C.y)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 3

# Coordinate system base vectors (C.i, C.j) are not being processed by ScalarField,
# so as a result of applying field they are kept untouched, like all other free variables in SymPy expression.
def test_dimensional_vector_to_field_conversion(test_args):
    field = sympy_vector_to_field(test_args.C.x * test_args.C.i + test_args.C.y * test_args.C.j)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == test_args.C.i + 2 * test_args.C.j

def test_empty_vector_to_field_conversion():
    field = sympy_vector_to_field(0)
    # applying empty field to a point results in zero value
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 0

def test_only_integer_vector_to_field_conversion():
    field = sympy_vector_to_field(1)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 1

def test_custom_names_vector_to_field_conversion():
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    field = sympy_vector_to_field(C1.r + 2 * C1.phi)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 5

def test_rotate_coordinates_vector_to_field_conversion(test_args):
    sympy_vector_field = test_args.C.x + test_args.C.y
    field = sympy_vector_to_field(sympy_vector_field)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 3
    theta = symbols("theta")
    B = test_args.C.orient_new_axis('B', theta, test_args.C.k)
    transformed_vector = express(sympy_vector_field, B, variables=True)
    result_transformed_field = sympy_vector_to_field(transformed_vector)
    assert transformed_vector == B.x * sin(theta) + B.x * cos(theta) - B.y * sin(theta) + B.y * cos(theta)
    assert result_transformed_field(field_point) == sin(theta) + cos(theta) - 2 * sin(theta) + 2 * cos(theta)

# when we express SymPy Vector to another coordinate system and base scalars (C.x, C.y) are
# not transformed, we have same base scalars as before.
def test_rotate_coordinates_without_variables_vector_to_field_conversion(test_args):
    theta = symbols("theta")
    B = test_args.C.orient_new_axis('B', theta, test_args.C.k)
    transformed_field = express(test_args.C.x + test_args.C.y, B)
    field = sympy_vector_to_field(transformed_field)
    field_point = FieldPoint(1, 2, 3)
    assert field(field_point) == 3

# Test apply_field()

# Result is a function that returns a scalar value at any point of the trajectory.
# Input field has X * Y, so as we are moving along X-axis or Y-axis of the trajectory,
# resulting value grows proportionally to X * Y.
def test_basic_field_apply(test_args):
    field = ScalarField(lambda p: p.y * p.x)
    # represents surface
    trajectory = [test_args.C.x, test_args.C.y]
    trajectory_value = apply_field(field, trajectory)
    assert trajectory_value == test_args.C.y * test_args.C.x

# Coordinate system is not necessary to apply field.
def test_parametrized_no_coord_system_field_apply():
    field = ScalarField(lambda p: -p.y)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_value = apply_field(field, trajectory)
    assert trajectory_value == -parameter

def test_parametrized_field_apply(test_args):
    result_field = sympy_vector_to_field(-test_args.C.y + 2 * test_args.C.x)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_value = apply_field(result_field, trajectory)
    assert trajectory_value == parameter

def test_sympy_field_apply(test_args):
    result_field = sympy_vector_to_field(-test_args.C.y + test_args.C.x)
    trajectory = [test_args.C.x, test_args.C.y]
    trajectory_value = apply_field(result_field, trajectory)
    assert trajectory_value == -test_args.C.y + test_args.C.x
