from collections import namedtuple
from pytest import fixture, raises
from sympy import atan, pi, sqrt, symbols, sin, cos
from sympy.vector import Vector as SympyVector, CoordSys3D, express
from symplyphysics.core.vectors.vectors import Vector, extract_coord_system_from_sympy_vector, sympy_vector_from_vector, vector_from_sympy_vector, vector_rebase


@fixture
def test_args():
    C = CoordSys3D("C")
    Args = namedtuple("Args", ["C"])
    return Args(C=C)


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

# Test Vector constructor

def test_basic_vector():
    vector = Vector([1, 2])
    assert vector.components == [1, 2]
    assert vector.coord_system == None

def test_coord_sys_vector(test_args):
    vector = Vector([1, 2], test_args.C)
    assert vector.components == [1, 2]
    assert vector.coord_system == test_args.C

def test_empty_vector():
    vector = Vector()
    assert vector.components == []
    assert vector.coord_system == None

# Test vector_from_sympy_vector()

def test_basic_sympy_to_array_conversion(test_args):
    vector = vector_from_sympy_vector(test_args.C.i + 2 * test_args.C.j)
    assert vector.components == [1, 2, 0]
    assert vector.coord_system == test_args.C

def test_order_sympy_to_array_conversion(test_args):
    vector = vector_from_sympy_vector(2 * test_args.C.j + test_args.C.i)
    assert vector.components == [1, 2, 0]

def test_skip_dimension_sympy_to_array_conversion(test_args):
    vector = vector_from_sympy_vector(test_args.C.i + 2 * test_args.C.k)
    assert vector.components == [1, 0, 2]

def test_empty_sympy_to_array_conversion():
    vector = vector_from_sympy_vector(SympyVector.zero)
    assert vector.components == []
    assert vector.coord_system == None

# If expression is not a SymPy Vector it will not be converted
def test_only_scalar_sympy_to_array_conversion(test_args):
    vector = vector_from_sympy_vector(test_args.C.x)
    assert vector.components == [test_args.C.x]
    assert vector.coord_system == None

def test_free_variable_sympy_to_array_conversion(test_args):
    x1 = symbols("x1")
    vector = vector_from_sympy_vector(test_args.C.i * x1)
    assert vector.components == [x1, 0, 0]
    assert vector.coord_system == test_args.C

def test_free_variable_empty_sympy_to_array_conversion():
    x1 = symbols("x1")
    vector = vector_from_sympy_vector(x1)
    assert vector.components == [x1]
    assert vector.coord_system == None

def test_custom_names_array_to_sympy_conversion():
    C1 = CoordSys3D("C1", vector_names=("r", "phi", "z"))
    vector = vector_from_sympy_vector(C1.r + 2 * C1.phi)
    assert vector.components == [1, 2, 0]
    assert vector.coord_system == C1

def test_rotate_coordinates_array_to_sympy_conversion(test_args):
    sympy_vector = test_args.C.i + test_args.C.j
    vector = vector_from_sympy_vector(sympy_vector)
    assert vector.components == [1, 1, 0]
    theta = symbols("theta")
    B = test_args.C.orient_new_axis('B', theta, test_args.C.k)
    sympy_transformed_vector = express(sympy_vector, B)
    assert sympy_transformed_vector == ((sin(theta) + cos(theta)) * B.i + (-sin(theta) + cos(theta)) * B.j)
    transformed_vector = vector_from_sympy_vector(sympy_transformed_vector)
    assert transformed_vector.components == [sin(theta) + cos(theta), -sin(theta) + cos(theta), 0]
    assert transformed_vector.coord_system == B

# Test sympy_vector_from_vector()

def test_custom_names_sympy_to_array_conversion():
    C1 = CoordSys3D("C1", vector_names=("r", "phi", "z"))
    sympy_vector = sympy_vector_from_vector(Vector([1, 2], C1))
    assert sympy_vector == C1.r + 2 * C1.phi

def test_multiple_coord_systems_sympy_to_array_conversion(test_args):
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    with raises(TypeError):
        vector_from_sympy_vector(test_args.C.i + 2 * C1.k)
    with raises(TypeError):
        vector_from_sympy_vector(test_args.C.i + C1.r * test_args.C.k)

def test_basic_array_to_sympy_conversion(test_args):
    sympy_vector = sympy_vector_from_vector(Vector([1, 2], test_args.C))
    assert sympy_vector == test_args.C.i + 2 * test_args.C.j

def test_skip_dimension_array_to_sympy_conversion(test_args):
    sympy_vector = sympy_vector_from_vector(Vector([1, 0, 2], test_args.C))
    assert sympy_vector == test_args.C.i + 2 * test_args.C.k

def test_4d_array_to_sympy_conversion(test_args):
    sympy_vector = sympy_vector_from_vector(Vector([1, 0, 2, 5], test_args.C))
    assert sympy_vector == test_args.C.i + 2 * test_args.C.k

def test_empty_array_to_sympy_conversion(test_args):
    sympy_vector = sympy_vector_from_vector(Vector([], test_args.C))
    assert sympy_vector == SympyVector.zero
    # only comparison with Vector.zero works
    assert sympy_vector != 0
    assert sympy_vector is not None

def test_empty_coord_sys_to_sympy_conversion():
    sympy_vector = sympy_vector_from_vector(Vector([1, 2]))
    assert sympy_vector == SympyVector.zero

def test_string_array_to_sympy_conversion(test_args):
    with raises(TypeError):
        sympy_vector_from_vector(Vector(["test"], test_args.C))

def test_rotate_coordinates_sympy_to_array_conversion(test_args):
    theta = symbols("theta")
    B = test_args.C.orient_new_axis('B', theta, test_args.C.k)
    sympy_vector = sympy_vector_from_vector(Vector([1, 2], B))
    assert sympy_vector == B.i + 2 * B.j
    transformed_vector = express(sympy_vector, test_args.C)
    assert transformed_vector == ((-2 * sin(theta) + cos(theta)) * test_args.C.i + (sin(theta) + 2 * cos(theta)) * test_args.C.j)

# Test vector_rebase()

def test_basic_vector_rebase(test_args):
    vector = Vector([test_args.C.x, test_args.C.y], test_args.C)

    # B is located at [1, 2] origin instead of [0, 0] of test_args.C
    B = test_args.C.locate_new('B', test_args.C.i + 2 * test_args.C.j)
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coord_system == B
    # Original field is not changed
    assert vector.coord_system == test_args.C
    assert vector_rebased.components == [B.x + 1, B.y + 2, 0]

# Simple numbers in vector are not scalars - they cannot be properly
# rebased to new coordinate system.
def test_plain_vector_rebase(test_args):
    vector = Vector([1, 2], test_args.C)

    # B is located at [1, 2] origin instead of [0, 0] of test_args.C
    B = test_args.C.locate_new('B', test_args.C.i + 2 * test_args.C.j)
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coord_system == B
    # Original field is not changed
    assert vector.coord_system == test_args.C
    assert vector_rebased.components == [1, 2, 0]

# Simple parameters in vector are not scalars - they cannot be properly
# rebased to new coordinate system.
def test_parameters_vector_rebase(test_args):
    parameter = symbols("parameter")
    vector = Vector([parameter, parameter], test_args.C)

    # B is located at [1, 2] origin instead of [0, 0] of test_args.C
    B = test_args.C.locate_new('B', test_args.C.i + 2 * test_args.C.j)
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coord_system == B
    # Original field is not changed
    assert vector.coord_system == test_args.C
    assert vector_rebased.components == [parameter, parameter, 0]

# Vector is not rebased if no original coordinate system was set.
def test_no_coord_system_vector_rebase(test_args):
    vector = Vector([test_args.C.x, test_args.C.y])
    assert vector.coord_system is None
    B = test_args.C.locate_new('B', test_args.C.i + 2 * test_args.C.j)
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coord_system == B
    assert vector_rebased.components == [test_args.C.x, test_args.C.y]

# Vector is not rebased if no target coordinate system was set.
def test_no_target_coord_system_vector_rebase(test_args):
    vector = Vector([test_args.C.x, test_args.C.y], test_args.C)
    assert vector.coord_system == test_args.C
    vector_rebased = vector_rebase(vector, None)
    assert vector_rebased.coord_system == None
    assert vector_rebased.components == [test_args.C.x, test_args.C.y]

# Rotation does not require vector defined with base scalars.
def test_rotate_vector_rebase(test_args):
    vector = Vector([1, 2], test_args.C)
    point = [1, 2]
    p1 = test_args.C.origin.locate_new('p1', point[0] * test_args.C.i + point[1] * test_args.C.j)
    p1_coordinates = p1.express_coordinates(test_args.C)
    assert p1_coordinates[0] == point[0]
    assert p1_coordinates[1] == point[1]

    B = test_args.C.orient_new_axis('B', pi/4, test_args.C.k)
    p1_coordinates_in_b = p1.express_coordinates(B)
    assert p1_coordinates_in_b[0] != point[0]

    transformed_point = [ p1_coordinates_in_b[0], p1_coordinates_in_b[1], 0 ]
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coord_system == B
    assert vector_rebased.components == [3 * sqrt(2) / 2, sqrt(2) / 2, 0]
    assert vector_rebased.components == transformed_point

# Test non-cartesian coordinate systems

#TODO: this test should pass if non-cartesian rebase is allowed. SymPy does not support such thing yet.
def test_spherical_vector_create(test_args):
    vector = Vector([1, 2], test_args.C)
    B = test_args.C.create_new('B', transformation="spherical")
    # vector should have r = sqrt(5) in polar coordinates
    # vector is in XY-plane, so theta angle should be pi/2
    # phi angle is atan(2/1)
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coord_system == B
    assert vector_rebased.components != [sqrt(5), pi/2, atan(2)]