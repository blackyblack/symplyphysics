from collections import namedtuple
from pytest import fixture, raises
from sympy import symbols
from sympy.vector import Vector, CoordSys3D
from symplyphysics.vectors import array_to_sympy_vector, sympy_vector_to_array


@fixture
def test_args():
    C = CoordSys3D('C')
    Args = namedtuple('Args', ['C'])
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

def test_free_variable_sympy_to_array_conversion(test_args):
    x1 = symbols('x1')
    result_array = sympy_vector_to_array(test_args.C.i * x1)
    assert result_array == [x1, 0, 0]

def test_free_variable_empty_sympy_to_array_conversion():
    x1 = symbols('x1')
    result_array = sympy_vector_to_array(x1)
    assert result_array == []

#TODO: add tests for CoordSys3D with custom vector names
#TODO: add tests for CoordSys3D with cylindrical and spherical coordinate systems
