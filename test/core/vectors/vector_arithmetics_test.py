from collections import namedtuple
from pytest import fixture, raises
from sympy import atan2, cos, pi, sin, sqrt
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem, coordinates_transform
from symplyphysics.core.vectors.vector_arithmetics import add_cartesian_vectors, cross_cartesian_vectors, dot_vectors, equal_vectors, scale_vector, vector_magnitude
from symplyphysics.core.vectors.vectors import Vector, vector_rebase
from test.test_decorators import unsupported_usage


@fixture(name="test_args")
def test_args_fixture():
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
    assert equal_vectors(Vector([1, 2, 3]), Vector([1, 2, 3]))
    assert equal_vectors(Vector([test_args.C.coord_system.x, 2]),
        Vector([test_args.C.coord_system.x, 2]))
    assert not equal_vectors(Vector([test_args.C.coord_system.x, 2]),
        Vector([test_args.C.coord_system.y, 2]))
    assert equal_vectors(Vector([test_args.C.coord_system.x - test_args.C.coord_system.x, 2]),
        Vector([0, 2]))
    assert equal_vectors(Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]))
    assert not equal_vectors(Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
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
    assert add_cartesian_vectors(Vector([1, 2]), Vector([1, 2])).components == [2, 4]
    assert add_cartesian_vectors(Vector([1, 2]), Vector([2, 1])).components == [3, 3]
    assert add_cartesian_vectors(Vector([1, 2 * 2]), Vector([1, 4])).components == [2, 8]
    assert add_cartesian_vectors(Vector([1, 0]), Vector([1])).components == [2, 0]
    assert add_cartesian_vectors(Vector([]), Vector([])).components == []
    assert add_cartesian_vectors(Vector([1, 2, 3]), Vector([1, 2, 3])).components == [2, 4, 6]
    assert add_cartesian_vectors(Vector([test_args.C.coord_system.x, 2]),
        Vector([test_args.C.coord_system.x, 2])).components == [2 * test_args.C.coord_system.x, 4]
    assert add_cartesian_vectors(Vector([test_args.C.coord_system.x, 2]),
        Vector([test_args.C.coord_system.y,
        2])).components == [test_args.C.coord_system.x + test_args.C.coord_system.y, 4]
    assert add_cartesian_vectors(
        Vector([test_args.C.coord_system.x - test_args.C.coord_system.x, 2]), Vector([1,
        2])).components == [1, 4]
    assert add_cartesian_vectors(
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.i,
        2])).components == [2 * test_args.C.coord_system.x * test_args.C.coord_system.i, 4]
    assert add_cartesian_vectors(
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
        Vector([test_args.C.coord_system.x * test_args.C.coord_system.j, 2])).components == [
        test_args.C.coord_system.x * test_args.C.coord_system.i +
        test_args.C.coord_system.x * test_args.C.coord_system.j, 4
        ]
    result_vector = add_cartesian_vectors(Vector([1, 2], test_args.C), Vector([1, 2], test_args.C))
    assert result_vector.components == [2, 4]
    assert result_vector.coordinate_system == test_args.C


def test_rebased_add_vectors(test_args):
    Cy = coordinates_transform(test_args.C, CoordinateSystem.System.CYLINDRICAL)
    # First vector has r = 5 and theta = 0
    vector1 = Vector([5, 0, 1], Cy)
    # Second vector has r = 5 and theta = pi/2
    vector2 = Vector([5, pi / 2, 1], Cy)
    # Their sum should be a vector with r = sqrt(50) and theta = pi/4
    vector1_cartesian = vector_rebase(vector1, test_args.C)
    vector2_cartesian = vector_rebase(vector2, test_args.C)
    vector_cartesian_sum = add_cartesian_vectors(vector1_cartesian, vector2_cartesian)
    vector_sum = vector_rebase(vector_cartesian_sum, Cy)
    assert vector_sum.components == [sqrt(50), pi / 4, 2]


def test_invalid_add_vectors(test_args):
    C1 = CoordinateSystem()
    with raises(TypeError):
        add_cartesian_vectors(Vector(["hello", 2]), Vector([2, 1]))
    with raises(TypeError):
        add_cartesian_vectors(Vector([1, 2]), Vector([[1, 2], 1]))
    with raises(AttributeError):
        add_cartesian_vectors(Vector([1, 2]), [1, 2])
    with raises(TypeError):
        add_cartesian_vectors(Vector([1, 2], test_args.C), Vector([1, 2]))
    with raises(TypeError):
        add_cartesian_vectors(Vector([1, 2], test_args.C), Vector([1, 2], C1))
    with raises(TypeError):
        add_cartesian_vectors(Vector([1, 2]), Vector([1, 2], test_args.C))

    # non-cartesian addition is not supported
    C2 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    with raises(ValueError):
        add_cartesian_vectors(Vector([1, 2], C2), Vector([1, 2], C2))
    C3 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)
    with raises(ValueError):
        add_cartesian_vectors(Vector([1, 2], C3), Vector([1, 2], C3))


@unsupported_usage
def test_strings_add_vectors():
    with raises(TypeError):
        add_cartesian_vectors(Vector(["hello", 2]), Vector([2, 1]))
    assert add_cartesian_vectors(Vector(["hello", 2]),
        Vector([" world"])).components == ["hello world", 2]


def test_basic_scale_vector(test_args):
    assert scale_vector(2, Vector([1, 2])).components == [2, 4]
    # Verify that magnitude has changed times scalar value, angle hasn't changed
    initial_vector_magnitude = sqrt(1**2 + 2**2)
    result_vector_magnitude = sqrt(2**2 + 4**2)
    assert result_vector_magnitude / initial_vector_magnitude == 2
    initial_vector_angle = atan2(2, 1)
    result_vector_angle = atan2(4, 2)
    assert initial_vector_angle == result_vector_angle

    assert scale_vector(2, Vector([2, 2])).components == [4, 4]
    # Additional check that magnitude has changed times scalar value
    initial_vector_magnitude = sqrt(2**2 + 2**2)
    result_vector_magnitude = sqrt(4**2 + 4**2)
    assert result_vector_magnitude / initial_vector_magnitude == 2

    # Try 3 dimensions
    assert scale_vector(2, Vector([1, 2, 3])).components == [2, 4, 6]
    # Verify that magnitude has changed times scalar value, angle hasn't changed
    initial_vector_magnitude = sqrt(1**2 + 2**2 + 3**2)
    result_vector_magnitude = sqrt(2**2 + 4**2 + 6**2)
    assert result_vector_magnitude / initial_vector_magnitude == 2
    initial_vector_angle = atan2(2, 1)
    result_vector_angle = atan2(4, 2)
    assert initial_vector_angle == result_vector_angle

    assert scale_vector(1, Vector([2, 1])).components == [2, 1]
    assert scale_vector(0.1, Vector([1, 4])).components == [0.1, 0.4]
    assert scale_vector(2, Vector([1])).components == [2]
    assert scale_vector(2, Vector([])).components == []
    assert scale_vector(0, Vector([1, 4])).components == [0, 0]
    assert scale_vector(-1, Vector([1, 4])).components == [-1, -4]
    assert scale_vector(2, Vector([test_args.C.coord_system.x,
        2])).components == [2 * test_args.C.coord_system.x, 4]
    assert scale_vector(test_args.C.coord_system.x, Vector([test_args.C.coord_system.y,
        2])).components == [
        test_args.C.coord_system.x * test_args.C.coord_system.y, test_args.C.coord_system.x * 2
        ]
    result_vector = scale_vector(2, Vector([1, 2], test_args.C))
    assert result_vector.components == [2, 4]
    assert result_vector.coordinate_system == test_args.C


def test_invalid_scale_vector():
    with raises(AttributeError):
        scale_vector(2, [1, 2])
    with raises(TypeError):
        scale_vector(Vector([1, 2]), Vector([1, 2]))


def test_cylindrical_scale_vector(test_args):
    vector1 = Vector([1, 2], test_args.C)
    Cy = coordinates_transform(test_args.C, CoordinateSystem.System.CYLINDRICAL)
    vector1_polar = vector_rebase(vector1, Cy)
    scaled_vector1_polar = scale_vector(5, vector1_polar)
    scaled_vector1 = vector_rebase(scaled_vector1_polar, test_args.C)

    # Verify that magnitude has changed times scalar value, angle hasn't changed
    initial_vector_magnitude = sqrt(vector1.components[0]**2 + vector1.components[1]**2)
    result_vector_magnitude = sqrt(scaled_vector1.components[0]**2 +
        scaled_vector1.components[1]**2)
    assert result_vector_magnitude / initial_vector_magnitude == 5
    initial_vector_angle = atan2(vector1.components[1], vector1.components[0])
    result_vector_angle = atan2(scaled_vector1.components[1], scaled_vector1.components[0])
    assert initial_vector_angle == result_vector_angle

    # Try 3 dimensions
    vector2 = Vector([1, 2, 3], test_args.C)
    Cy = coordinates_transform(test_args.C, CoordinateSystem.System.CYLINDRICAL)
    vector2_polar = vector_rebase(vector2, Cy)
    scaled_vector2_polar = scale_vector(5, vector2_polar)
    scaled_vector2 = vector_rebase(scaled_vector2_polar, test_args.C)

    # Verify that magnitude has changed times scalar value, angle hasn't changed
    initial_vector2_magnitude = sqrt(vector2.components[0]**2 + vector2.components[1]**2 +
        vector2.components[2]**2)
    result_vector2_magnitude = sqrt(scaled_vector2.components[0]**2 +
        scaled_vector2.components[1]**2 + scaled_vector2.components[2]**2)
    assert result_vector2_magnitude / initial_vector2_magnitude == 5
    # It's hard to verify angle so only validate magnitude

    assert scale_vector(5, Vector([1, 2], Cy)).components == [5, 2]
    assert scale_vector(5, Vector([1, 2, 3], Cy)).components == [5, 2, 15]


def test_spherical_scale_vector(test_args):
    vector1 = Vector([1, 2], test_args.C)
    Cs = coordinates_transform(test_args.C, CoordinateSystem.System.SPHERICAL)
    vector1_polar = vector_rebase(vector1, Cs)
    scaled_vector1_polar = scale_vector(5, vector1_polar)
    scaled_vector1 = vector_rebase(scaled_vector1_polar, test_args.C)

    # Verify that magnitude has changed times scalar value, angle hasn't changed
    initial_vector_magnitude = sqrt(vector1.components[0]**2 + vector1.components[1]**2)
    result_vector_magnitude = sqrt(scaled_vector1.components[0]**2 +
        scaled_vector1.components[1]**2)
    assert result_vector_magnitude / initial_vector_magnitude == 5
    initial_vector_angle = atan2(vector1.components[1], vector1.components[0])
    result_vector_angle = atan2(scaled_vector1.components[1], scaled_vector1.components[0])
    assert initial_vector_angle == result_vector_angle

    # Try 3 dimensions
    vector2 = Vector([1, 2, 3], test_args.C)
    Cs = coordinates_transform(test_args.C, CoordinateSystem.System.SPHERICAL)
    vector2_polar = vector_rebase(vector2, Cs)
    scaled_vector2_polar = scale_vector(5, vector2_polar)
    scaled_vector2 = vector_rebase(scaled_vector2_polar, test_args.C)

    # Verify that magnitude has changed times scalar value, angle hasn't changed
    initial_vector2_magnitude = sqrt(vector2.components[0]**2 + vector2.components[1]**2 +
        vector2.components[2]**2)
    result_vector2_magnitude = sqrt(scaled_vector2.components[0]**2 +
        scaled_vector2.components[1]**2 + scaled_vector2.components[2]**2)
    assert result_vector2_magnitude / initial_vector2_magnitude == 5
    # It's hard to verify angle so only validate magnitude

    assert scale_vector(5, Vector([1, 2], Cs)).components == [5, 2]
    assert scale_vector(5, Vector([1, 2, 3], Cs)).components == [5, 2, 3]


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
    assert dot_vectors(Vector([1, 2, 3]), Vector([3, 2, 1])) == 10
    assert dot_vectors(Vector([test_args.C.coord_system.x, 2]),
        Vector([test_args.C.coord_system.x, 2])) == test_args.C.coord_system.x**2 + 4
    assert dot_vectors(Vector([test_args.C.coord_system.x, 2]),
        Vector([test_args.C.coord_system.y,
        2])) == test_args.C.coord_system.x * test_args.C.coord_system.y + 4
    assert dot_vectors(Vector([test_args.C.coord_system.x - test_args.C.coord_system.x, 2]),
        Vector([1, 2])) == 4
    assert dot_vectors(Vector([1, 2], test_args.C), Vector([1, 2], test_args.C)) == 5
    # Verify that orthogonal vectors have zero dot product
    assert dot_vectors(Vector([1, 0]), Vector([0, 1])) == 0
    initial_vector_angle = atan2(3, 2)
    result_vector_angle = initial_vector_angle + pi / 2
    orthogonal_unit_vector = Vector([cos(result_vector_angle), sin(result_vector_angle)])
    assert dot_vectors(Vector([2, 3]), orthogonal_unit_vector) == 0
    # Verify that parallel vectors have dot product equal to their magnitudes multiplication
    first_vector_magnitude = sqrt(1**2 + 2**2)
    second_vector_magnitude = sqrt(2**2 + 4**2)
    assert dot_vectors(Vector([1, 2]), Vector([2,
        4])) == first_vector_magnitude * second_vector_magnitude


def test_invalid_dot_product(test_args):
    with raises(ValueError):
        dot_vectors(Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
            Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]))
    with raises(ValueError):
        dot_vectors(Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
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


def test_non_cartesian_dot_product():
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    assert dot_vectors(Vector([1, 0], C1), Vector([1, pi / 2], C1)) == 0
    assert dot_vectors(Vector([1, 0], C1), Vector([1, pi / 4], C1)) == (sqrt(2) / 2)
    assert dot_vectors(Vector([1, 0, 0], C1), Vector([1, pi / 2, 0], C1)) == 0
    assert dot_vectors(Vector([1, 0, 0], C1), Vector([0, 0, 1], C1)) == 0

    # When r = 0 and z > 0, vector is directed up and is immune to rotations,
    # therefore rotated vector is parallel to original
    assert dot_vectors(Vector([0, 0, 1], C1), Vector([0, pi / 2, 1], C1)) == 1

    C2 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)
    assert dot_vectors(Vector([1, 0], C2), Vector([1, pi / 2], C2)) == 0
    assert dot_vectors(Vector([1, 0], C2), Vector([1, pi / 4], C2)) == (sqrt(2) / 2)
    assert dot_vectors(Vector([1, 0, 0], C2), Vector([1, pi / 2, 0], C2)) == 0
    assert dot_vectors(Vector([1, pi / 2, 0], C2), Vector([1, pi / 2, pi / 2], C2)) == 0

    # When theta angle = 0, vector is directed up and is immune to rotations,
    # therefore rotated vector is parallel to original
    assert dot_vectors(Vector([1, 0, 0], C2), Vector([1, 0, pi / 2], C2)) == 1


def test_basic_magnitude_vectors(test_args):
    assert vector_magnitude(Vector([1, 2])) == sqrt(5)
    assert vector_magnitude(Vector([])) == 0
    assert vector_magnitude(Vector([1, 2, 3])) == sqrt(14)
    assert vector_magnitude(Vector([test_args.C.coord_system.x,
        2])) == sqrt(test_args.C.coord_system.x**2 + 4)


def test_non_cartesian_magnitude_vectors():
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    assert vector_magnitude(Vector([3, 0], C1)) == 3
    assert vector_magnitude(Vector([3, 0, 2], C1)) == sqrt(13)
    assert vector_magnitude(Vector([3, pi / 4, 2], C1)) == sqrt(13)

    C2 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)
    assert vector_magnitude(Vector([3, 0], C2)) == 3
    assert vector_magnitude(Vector([3, 0, 2], C2)) == 3
    assert vector_magnitude(Vector([3, pi / 4, 2], C2)) == 3


def test_basic_cross_product(test_args):
    # Parallel vectors have zero cross product
    assert cross_cartesian_vectors(Vector([1, 2]), Vector([1, 2])).components == [0, 0, 0]
    assert cross_cartesian_vectors(Vector([1, 2]), Vector([2, 1])).components == [0, 0, -3]
    assert cross_cartesian_vectors(Vector([1, 2 * 2]), Vector([1, 4])).components == [0, 0, 0]
    assert cross_cartesian_vectors(Vector([1, 0]), Vector([1])).components == [0, 0, 0]
    assert cross_cartesian_vectors(Vector([]), Vector([])).components == [0, 0, 0]
    assert cross_cartesian_vectors(Vector([1, 2, 3]), Vector([3, 2, 1])).components == [-4, 8, -4]
    assert cross_cartesian_vectors(Vector([test_args.C.coord_system.x, 2]),
        Vector([test_args.C.coord_system.x, 2])).components == [0, 0, 0]
    assert cross_cartesian_vectors(Vector([test_args.C.coord_system.x, 2]),
        Vector([test_args.C.coord_system.y,
        2])).components == [0, 0, 2 * test_args.C.coord_system.x - 2 * test_args.C.coord_system.y]
    assert cross_cartesian_vectors(
        Vector([test_args.C.coord_system.x - test_args.C.coord_system.x, 2]), Vector([1,
        2])).components == [0, 0, -2]
    assert cross_cartesian_vectors(Vector([1, 2], test_args.C), Vector([1, 2],
        test_args.C)).components == [0, 0, 0]
    # Verify that cross product vector is orthogonal to both input vectors
    input_vector1 = Vector([1, 0])
    input_vector2 = Vector([0, 1])
    result_cross = cross_cartesian_vectors(input_vector1, input_vector2)
    assert result_cross.components == [0, 0, 1]
    assert dot_vectors(result_cross, input_vector1) == 0
    assert dot_vectors(result_cross, input_vector2) == 0

    # For orthogonal vectors cross product magnitude equals to input vector magnitudes multiplication
    input_vector1 = Vector([2, 0])
    input_vector2 = Vector([0, 3])
    result_cross = cross_cartesian_vectors(input_vector1, input_vector2)
    assert vector_magnitude(
        result_cross) == vector_magnitude(input_vector1) * vector_magnitude(input_vector2)

    # For cross product magnitude equals to area of parallelogram spanned by input vectors
    input_vector1 = Vector([2, 3, 4])
    input_vector2 = Vector([3, 4, 5])
    result_cross = cross_cartesian_vectors(input_vector1, input_vector2)
    assert vector_magnitude(
        result_cross
    )**2 == vector_magnitude(input_vector1)**2 * vector_magnitude(input_vector2)**2 - dot_vectors(
        input_vector1, input_vector2)**2


def test_invalid_cross_product(test_args):
    with raises(TypeError):
        cross_cartesian_vectors(
            Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
            Vector([2, test_args.C.coord_system.x * test_args.C.coord_system.i]))
    with raises(TypeError):
        cross_cartesian_vectors(
            Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
            Vector([2, test_args.C.coord_system.x * test_args.C.coord_system.j]))
    C1 = CoordinateSystem()
    with raises(TypeError):
        cross_cartesian_vectors(Vector(["hello", 2]), Vector([2, 1]))
    with raises(TypeError):
        cross_cartesian_vectors(Vector([2, 2]), Vector(["hello", "world"]))
    with raises(TypeError):
        cross_cartesian_vectors(Vector([1, 2]), Vector([[1, 2], 1]))
    with raises(AttributeError):
        cross_cartesian_vectors(Vector([1, 2]), [1, 2])
    with raises(TypeError):
        cross_cartesian_vectors(Vector([1, 2], test_args.C), Vector([1, 2]))
    with raises(TypeError):
        cross_cartesian_vectors(Vector([1, 2], test_args.C), Vector([1, 2], C1))
    with raises(TypeError):
        cross_cartesian_vectors(Vector([1, 2]), Vector([1, 2], test_args.C))
    with raises(ValueError):
        cross_cartesian_vectors(Vector([1, 2, 3, 4], test_args.C), Vector([1, 2], test_args.C))

    # non-cartesian addition is not supported
    C2 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    with raises(ValueError):
        cross_cartesian_vectors(Vector([1, 2], C2), Vector([1, 2], C2))
    C3 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)
    with raises(ValueError):
        cross_cartesian_vectors(Vector([1, 2], C3), Vector([1, 2], C3))


def test_rebased_cross_product(test_args):
    Cy = coordinates_transform(test_args.C, CoordinateSystem.System.CYLINDRICAL)
    # First vector has r = 5 and theta = 0
    vector1 = Vector([5, 0, 1], Cy)
    # Second vector has r = 5 and theta = pi/2
    vector2 = Vector([5, pi / 2, 1], Cy)
    # Their product should be a vector with r = sqrt(50), theta = -3 * pi/4, z = 25
    vector1_cartesian = vector_rebase(vector1, test_args.C)
    vector2_cartesian = vector_rebase(vector2, test_args.C)
    vector_cartesian_cross = cross_cartesian_vectors(vector1_cartesian, vector2_cartesian)
    vector_cross = vector_rebase(vector_cartesian_cross, Cy)
    assert vector_cross.components == [sqrt(50), -3 * pi / 4, 25]
