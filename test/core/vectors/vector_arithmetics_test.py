from collections import namedtuple
from pytest import fixture, raises
from sympy import atan, cos, pi, sin, sqrt
from sympy.vector import  CoordSys3D
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.vectors.vector_arithmetics import add_vectors, dot_vectors, equal_vectors, scale_vector
from symplyphysics.core.vectors.vectors import Vector
from test.test_decorators import unsupported_usage


@fixture
def test_args():
    C = CoordinateSystem()
    Args = namedtuple("Args", ["C"])
    return Args(C=C)


def test_basic_equal_vectors(test_args):
    assert equal_vectors(Vector([1, 2]), Vector([1, 2]))
    assert not equal_vectors(Vector([1, 2]), Vector([2, 1]))
    assert equal_vectors(Vector([1, 2 * 2]), Vector([1, 4]))
    assert equal_vectors(Vector([1, 0]), Vector([1]))
    assert equal_vectors(Vector([]), Vector([]))
    assert equal_vectors(Vector([0, 0]), Vector([]))
    assert equal_vectors(Vector([]), Vector([0, 0]))
    assert equal_vectors(Vector([test_args.C.coord_system.x, 2]), Vector([test_args.C.coord_system.x, 2]))
    assert not equal_vectors(Vector([test_args.C.coord_system.x, 2]), Vector([test_args.C.coord_system.y, 2]))
    assert equal_vectors(Vector([test_args.C.coord_system.x - test_args.C.coord_system.x, 2]), Vector([0, 2]))
    assert equal_vectors(
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]))
    assert not equal_vectors(
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.j, 2]))
    assert equal_vectors(Vector([1, 2], test_args.C), Vector([1, 2], test_args.C))

def test_invalid_equal_vectors(test_args):
    C1 = CoordinateSystem()
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
    assert add_vectors(Vector([test_args.C.coord_system.x, 2]), Vector([test_args.C.coord_system.x, 2])).components == [2 * test_args.C.coord_system.x, 4]
    assert add_vectors(Vector([test_args.C.coord_system.x, 2]), Vector([test_args.C.coord_system.y, 2])).components == [test_args.C.coord_system.x + test_args.C.coord_system.y, 4]
    assert add_vectors(
        Vector([test_args.C.coord_system.x - test_args.C.coord_system.x, 2]),
        Vector([1, 2])).components == [1, 4]
    assert add_vectors(
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2])
        ).components == [2 * test_args.C.coord_system.x * test_args.C.coord_system.i, 4]
    assert add_vectors(
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.j, 2])
        ).components == [test_args.C.coord_system.x * test_args.C.coord_system.i + test_args.C.coord_system.x * test_args.C.coord_system.j, 4]
    result_vector = add_vectors(Vector([1, 2], test_args.C), Vector([1, 2], test_args.C))
    assert result_vector.components == [2, 4]
    assert result_vector.coordinate_system == test_args.C

def test_invalid_add_vectors(test_args):
    C1 = CoordinateSystem()
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

def test_basic_scale_vector(test_args):
    assert scale_vector(2, Vector([1, 2])).components == [2, 4]
    # Verify that magnitude has changed times scalar value, angle hasn't changed
    initial_vector_magnitude = sqrt(1**2 + 2**2)
    result_vector_magnitude = sqrt(2**2 + 4**2)
    assert result_vector_magnitude / initial_vector_magnitude == 2
    initial_vector_angle = atan(2 / 1)
    result_vector_angle = atan(4 / 2)
    assert initial_vector_angle == result_vector_angle
    assert scale_vector(1, Vector([2, 1])).components == [2, 1]
    assert scale_vector(0.1, Vector([1, 4])).components == [0.1, 0.4]
    assert scale_vector(2, Vector([1])).components == [2]
    assert scale_vector(2, Vector([])).components == []
    assert scale_vector(0, Vector([1, 4])).components == [0, 0]
    assert scale_vector(-1, Vector([1, 4])).components == [-1, -4]
    assert scale_vector(2, Vector([test_args.C.coord_system.x, 2])).components == [2 * test_args.C.coord_system.x, 4]
    assert scale_vector(
        test_args.C.coord_system.x,
        Vector([test_args.C.coord_system.y, 2])
        ).components == [test_args.C.coord_system.x * test_args.C.coord_system.y, test_args.C.coord_system.x * 2]
    result_vector = scale_vector(2, Vector([1, 2], test_args.C))
    assert result_vector.components == [2, 4]
    assert result_vector.coordinate_system == test_args.C

def test_invalid_scale_vector():
    with raises(AttributeError):
        scale_vector(2, [1, 2])
    with raises(TypeError):
        scale_vector(Vector([1, 2]), Vector([1, 2]))

@unsupported_usage
def test_unsupported_scale_vector():
    assert scale_vector(2, Vector(["hello", 2])).components == ["hellohello", 4]
    assert scale_vector(2, Vector([[1, 2], 1])).components == [[1, 2, 1, 2], 2]
    assert scale_vector([1, 2], Vector([1, 2])).components == [[1, 2], [1, 2, 1, 2]]

def test_basic_dot_product(test_args):
    assert dot_vectors(Vector([1, 2]), Vector([1, 2])) == 5
    assert dot_vectors(Vector([1, 2]), Vector([2, 1])) == 4
    assert dot_vectors(Vector([1, 2 * 2]), Vector([1, 4])) == 17
    assert dot_vectors(Vector([1, 0]), Vector([1])) == 1
    assert dot_vectors(Vector([]), Vector([])) == 0
    assert dot_vectors(
        Vector([test_args.C.coord_system.x, 2]),
        Vector([test_args.C.coord_system.x, 2])
        ) == test_args.C.coord_system.x**2 + 4
    assert dot_vectors(
        Vector([test_args.C.coord_system.x, 2]),
        Vector([test_args.C.coord_system.y, 2])
        ) == test_args.C.coord_system.x * test_args.C.coord_system.y + 4
    assert dot_vectors(
        Vector([test_args.C.coord_system.x - test_args.C.coord_system.x, 2]),
        Vector([1, 2])
        ) == 4
    assert dot_vectors(
        Vector([1, 2], test_args.C),
        Vector([1, 2], test_args.C)
        ) == 5
    # Verify that orthogonal vectors have zero dot product
    assert dot_vectors(Vector([1, 0]), Vector([0, 1])) == 0
    initial_vector_angle = atan(3 / 2)
    result_vector_angle = initial_vector_angle + pi / 2
    orthogonal_unit_vector = Vector([cos(result_vector_angle), sin(result_vector_angle)])
    assert dot_vectors(Vector([2, 3]), orthogonal_unit_vector) == 0
    # Verify that parallel vectors have dot product equal to their magnitudes multiplication
    first_vector_magnitude = sqrt(1**2 + 2**2)
    second_vector_magnitude = sqrt(2**2 + 4**2)
    assert dot_vectors(Vector([1, 2]), Vector([2, 4])) == first_vector_magnitude * second_vector_magnitude

def test_invalid_dot_product(test_args):
    with raises(ValueError):
        dot_vectors(
            Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
            Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]))
    with raises(ValueError):
        dot_vectors(
            Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
            Vector([test_args.C.coord_system.x * test_args.C.coord_system.j, 2]))
    C1 = CoordinateSystem()
    with raises(TypeError):
        dot_vectors(Vector(["hello", 2]), Vector([2, 1]))
    with raises(TypeError):
        dot_vectors(Vector([2, 2]), Vector(["hello", "world"]))
    with raises(TypeError):
        dot_vectors(Vector([1, 2]), Vector([[1, 2], 1]))
    with raises(AttributeError):
        dot_vectors(Vector([1, 2]), [1, 2])
    with raises(TypeError):
        dot_vectors(Vector([1, 2], test_args.C), Vector([1, 2]))
    with raises(TypeError):
        dot_vectors(Vector([1, 2], test_args.C), Vector([1, 2], C1))
    with raises(TypeError):
        dot_vectors(Vector([1, 2]), Vector([1, 2], test_args.C))

#TODO: non-cartesian coordinate systems are not supported. For example, scalar multiplication
#      of vector in polar coordinates [r, phi] should change it's magnitude but should not
#      change it's direction. So only r should be mutiplied and phi should stay the same.
#TODO: same for Dot product and vector length formulas.
#TODO: possible solution is to convert to SymPy vectors, make calculations and convert back to
#      Vector.
