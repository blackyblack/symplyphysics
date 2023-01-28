from collections import namedtuple
from pytest import fixture, raises
from sympy.vector import  CoordSys3D
from symplyphysics.core.vectors.vector_arithmetics import add_vectors, equal_vectors
from symplyphysics.core.vectors.vectors import Vector
from test.test_decorators import unsupported_usage


@fixture
def test_args():
    C = CoordSys3D("C")
    Args = namedtuple("Args", ["C"])
    return Args(C=C)


def test_basic_equal_vectors(test_args):
    assert equal_vectors(Vector([1, 2]), Vector([1, 2]))
    assert not equal_vectors(Vector([1, 2]), Vector([2, 1]))
    assert equal_vectors(Vector([1, 2 * 2]), Vector([1, 4]))
    assert equal_vectors(Vector([1, 0]), Vector([1]))
    assert equal_vectors(Vector([]), Vector([]))
    assert equal_vectors(Vector([test_args.C.x, 2]), Vector([test_args.C.x, 2]))
    assert not equal_vectors(Vector([test_args.C.x, 2]), Vector([test_args.C.y, 2]))
    assert equal_vectors(Vector([test_args.C.x - test_args.C.x, 2]), Vector([0, 2]))
    assert equal_vectors(Vector([test_args.C.x * test_args.C.i, 2]), Vector([test_args.C.x * test_args.C.i, 2]))
    assert not equal_vectors(Vector([test_args.C.x * test_args.C.i, 2]), Vector([test_args.C.x * test_args.C.j, 2]))
    assert equal_vectors(Vector([1, 2], test_args.C), Vector([1, 2], test_args.C))

def test_invalid_equal_vectors(test_args):
    C1 = CoordSys3D("C1")
    with raises(TypeError):
        equal_vectors(Vector(["hello", 2]), Vector([2, 1]))
    with raises(AttributeError):
        equal_vectors(Vector([1, 2]), [1, 2])
    with raises(TypeError):
        equal_vectors(Vector([1, 2], test_args.C), Vector([1, 2]))
    with raises(TypeError):
        equal_vectors(Vector([1, 2], test_args.C), Vector([1, 2], C1))
    with raises(TypeError):
        equal_vectors(Vector([1, 2]), Vector([1, 2], test_args.C))

def test_basic_add_vectors(test_args):
    assert add_vectors(Vector([1, 2]), Vector([1, 2])).components == [2, 4]
    assert add_vectors(Vector([1, 2]), Vector([2, 1])).components == [3, 3]
    assert add_vectors(Vector([1, 2 * 2]), Vector([1, 4])).components == [2, 8]
    assert add_vectors(Vector([1, 0]), Vector([1])).components == [2, 0]
    assert add_vectors(Vector([]), Vector([])).components == []
    assert add_vectors(Vector([test_args.C.x, 2]), Vector([test_args.C.x, 2])).components == [2* test_args.C.x, 4]
    assert add_vectors(Vector([test_args.C.x, 2]), Vector([test_args.C.y, 2])).components == [test_args.C.x + test_args.C.y, 4]
    assert add_vectors(Vector([test_args.C.x - test_args.C.x, 2]), Vector([0, 2])).components == [0, 4]
    assert add_vectors(Vector([test_args.C.x * test_args.C.i, 2]), Vector([test_args.C.x * test_args.C.i, 2])).components == [2 * test_args.C.x * test_args.C.i, 4]
    assert add_vectors(Vector([test_args.C.x * test_args.C.i, 2]), Vector([test_args.C.x * test_args.C.j, 2])).components == [test_args.C.x * test_args.C.i + test_args.C.x * test_args.C.j, 4]
    result_vector = add_vectors(Vector([1, 2], test_args.C), Vector([1, 2], test_args.C))
    assert result_vector.components == [2, 4]
    assert result_vector.coord_system == test_args.C

def test_invalid_add_vectors(test_args):
    C1 = CoordSys3D("C1")
    with raises(TypeError):
        add_vectors(Vector(["hello", 2]), Vector([2, 1]))
    with raises(TypeError):
        add_vectors(Vector([1, 2]), Vector([[1, 2], 1]))
    with raises(AttributeError):
        add_vectors(Vector([1, 2]), [1, 2])
    with raises(TypeError):
        add_vectors(Vector([1, 2], test_args.C), Vector([1, 2]))
    with raises(TypeError):
        add_vectors(Vector([1, 2], test_args.C), Vector([1, 2], C1))
    with raises(TypeError):
        add_vectors(Vector([1, 2]), Vector([1, 2], test_args.C))

@unsupported_usage
def test_strings_add_vectors():
    with raises(TypeError):
        add_vectors(Vector(["hello", 2]), Vector([2, 1]))
    assert add_vectors(Vector(["hello", 2]), Vector([" world"])).components == ["hello world", 2]

