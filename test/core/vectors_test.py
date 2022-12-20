from collections import namedtuple
from pytest import fixture, raises
from sympy import symbols, sin, cos
from sympy.vector import Vector, CoordSys3D, express
from symplyphysics.core.vectors import array_to_sympy_vector, sympy_vector_to_array


@fixture
def test_args():
    C = CoordSys3D("C")
    Args = namedtuple("Args", ["C"])
    return Args(C=C)


def test_basic_array_to_sympy_conversion(test_args):
    sympy_vector = array_to_sympy_vector(test_args.C, [1, 2])
    assert sympy_vector == test_args.C.i + 2 * test_args.C.j

def test_skip_dimension_array_to_sympy_conversion(test_args):
    sympy_vector = array_to_sympy_vector(test_args.C, [1, 0, 2])
    assert sympy_vector == test_args.C.i + 2 * test_args.C.k

def test_4d_array_to_sympy_conversion(test_args):
    sympy_vector = array_to_sympy_vector(test_args.C, [1, 0, 2, 5])
    assert sympy_vector == test_args.C.i + 2 * test_args.C.k

def test_empty_array_to_sympy_conversion(test_args):
    sympy_vector = array_to_sympy_vector(test_args.C, [])
    assert sympy_vector == Vector.zero
    # only comparison with Vector.zero works
    assert sympy_vector != 0
    assert sympy_vector is not None

def test_string_array_to_sympy_conversion(test_args):
    with raises(TypeError):
        array_to_sympy_vector(test_args.C, ["test"])

def test_basic_sympy_to_array_conversion(test_args):
    result_array = sympy_vector_to_array(test_args.C.i + 2 * test_args.C.j)
    assert result_array == [1, 2, 0]

def test_order_sympy_to_array_conversion(test_args):
    result_array = sympy_vector_to_array(2 * test_args.C.j + test_args.C.i)
    assert result_array == [1, 2, 0]

def test_skip_dimension_sympy_to_array_conversion(test_args):
    result_array = sympy_vector_to_array(test_args.C.i + 2 * test_args.C.k)
    assert result_array == [1, 0, 2]

def test_empty_sympy_to_array_conversion():
    result_array = sympy_vector_to_array(Vector.zero)
    assert result_array == []

# Expression should represent vector, so it should contain one of coordinate system base vectors, eg C.i
def test_only_scalar_sympy_to_array_conversion(test_args):
    with raises(TypeError):
        sympy_vector_to_array(test_args.C.x)

def test_free_variable_sympy_to_array_conversion(test_args):
    x1 = symbols("x1")
    result_array = sympy_vector_to_array(test_args.C.i * x1)
    assert result_array == [x1, 0, 0]

def test_free_variable_empty_sympy_to_array_conversion():
    x1 = symbols("x1")
    result_array = sympy_vector_to_array(x1)
    assert result_array == []

def test_custom_names_array_to_sympy_conversion():
    C1 = CoordSys3D("C1", vector_names=("r", "phi", "z"))
    result_array = sympy_vector_to_array(C1.r + 2 * C1.phi)
    assert result_array == [1, 2, 0]

def test_custom_names_sympy_to_array_conversion():
    C1 = CoordSys3D("C1", vector_names=("r", "phi", "z"))
    sympy_vector = array_to_sympy_vector(C1, [1, 2])
    assert sympy_vector == C1.r + 2 * C1.phi

def test_rotate_coordinates_array_to_sympy_conversion(test_args):
    sympy_vector = test_args.C.i + test_args.C.j
    result_array = sympy_vector_to_array(sympy_vector)
    assert result_array == [1, 1, 0]
    theta = symbols("theta")
    B = test_args.C.orient_new_axis('B', theta, test_args.C.k)
    transformed_vector = express(sympy_vector, B)
    assert transformed_vector == ((sin(theta) + cos(theta)) * B.i + (-sin(theta) + cos(theta)) * B.j)
    result_transformed_array = sympy_vector_to_array(transformed_vector)
    assert result_transformed_array == [sin(theta) + cos(theta), -sin(theta) + cos(theta), 0]

def test_rotate_coordinates_sympy_to_array_conversion(test_args):
    theta = symbols("theta")
    B = test_args.C.orient_new_axis('B', theta, test_args.C.k)
    sympy_vector = array_to_sympy_vector(B, [1, 2])
    assert sympy_vector == B.i + 2 * B.j
    transformed_vector = express(sympy_vector, test_args.C)
    assert transformed_vector == ((-2 * sin(theta) + cos(theta)) * test_args.C.i + (sin(theta) + 2 * cos(theta)) * test_args.C.j)

def test_multiple_coord_systems_sympy_to_array_conversion(test_args):
    C1 = CoordSys3D("C1", variable_names=("r", "phi", "z"))
    with raises(TypeError):
        sympy_vector_to_array(test_args.C.i + 2 * C1.k)
    with raises(TypeError):
        sympy_vector_to_array(test_args.C.i + C1.r * test_args.C.k)
