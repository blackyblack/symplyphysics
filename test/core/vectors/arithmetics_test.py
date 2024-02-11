from collections import namedtuple
from pytest import fixture, raises
from sympy import atan2, cos, pi, sin, sqrt
from symplyphysics import (SI, Quantity, dimensionless, units, QuantityVector, Vector,
    CoordinateSystem, coordinates_transform)
from symplyphysics.core.vectors.arithmetics import (add_cartesian_quantity_vectors,
    add_cartesian_vectors, cross_cartesian_quantity_vectors, cross_cartesian_vectors, dot_vectors,
    equal_vectors, quantity_vector_magnitude, quantity_vector_unit, scale_quantity_vector,
    scale_vector, vector_magnitude, dot_quantity_vectors, vector_unit)

# pylint: disable=too-many-locals

Args = namedtuple("Args", ["C"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C = CoordinateSystem()
    return Args(C=C)


def test_basic_equal_vectors(test_args: Args) -> None:
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


def test_invalid_equal_vectors(test_args: Args) -> None:
    C1 = CoordinateSystem()
    with raises(TypeError):
        equal_vectors(Vector([1, 2], test_args.C), Vector([1, 2], C1))


def test_basic_add_vectors(test_args: Args) -> None:
    assert add_cartesian_vectors(Vector([1, 2]), Vector([1, 2])).components == [2, 4]
    assert add_cartesian_vectors(Vector([1, 2]), Vector([2, 1])).components == [3, 3]
    assert add_cartesian_vectors(Vector([1, 2 * 2]), Vector([1, 4])).components == [2, 8]
    assert add_cartesian_vectors(Vector([1, 0]), Vector([1])).components == [2, 0]
    assert len(add_cartesian_vectors(Vector([]), Vector([])).components) == 0
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


def test_rebased_add_vectors(test_args: Args) -> None:
    Cy = coordinates_transform(test_args.C, CoordinateSystem.System.CYLINDRICAL)
    # First vector has r = 5 and theta = 0
    vector1 = Vector([5, 0, 1], Cy)
    # Second vector has r = 5 and theta = pi/2
    vector2 = Vector([5, pi / 2, 1], Cy)
    # Their sum should be a vector with r = sqrt(50) and theta = pi/4
    vector1_cartesian = vector1.rebase(test_args.C)
    vector2_cartesian = vector2.rebase(test_args.C)
    vector_cartesian_sum = add_cartesian_vectors(vector1_cartesian, vector2_cartesian)
    vector_sum = vector_cartesian_sum.rebase(Cy)
    assert vector_sum.components == [sqrt(50), pi / 4, 2]


def test_invalid_add_vectors(test_args: Args) -> None:
    C1 = CoordinateSystem()
    with raises(TypeError):
        add_cartesian_vectors(Vector([1, 2], test_args.C), Vector([1, 2], C1))
    # non-cartesian addition is not supported
    C2 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    with raises(ValueError):
        add_cartesian_vectors(Vector([1, 2], C2), Vector([1, 2], C2))
    C3 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)
    with raises(ValueError):
        add_cartesian_vectors(Vector([1, 2], C3), Vector([1, 2], C3))


def test_basic_scale_vector(test_args: Args) -> None:
    initial = Vector([1, 2])
    result = scale_vector(2, initial).components
    assert result == [2, 4]
    # Verify that magnitude has changed times scalar value, angle hasn't changed
    initial_vector_magnitude = sqrt(initial.components[0]**2 + initial.components[1]**2)
    result_vector_magnitude = sqrt(result[0]**2 + result[1]**2)
    assert result_vector_magnitude / initial_vector_magnitude == 2
    initial_vector_angle = atan2(initial.components[1], initial.components[0])
    result_vector_angle = atan2(result[1], result[0])
    assert initial_vector_angle == result_vector_angle

    initial = Vector([2, 2])
    result = scale_vector(2, initial).components
    assert result == [4, 4]
    # Additional check that magnitude has changed times scalar value
    initial_vector_magnitude = sqrt(initial.components[0]**2 + initial.components[1]**2)
    result_vector_magnitude = sqrt(result[0]**2 + result[1]**2)
    assert result_vector_magnitude / initial_vector_magnitude == 2

    # Try 3 dimensions
    initial = Vector([1, 2, 3])
    result = scale_vector(2, initial).components
    assert result == [2, 4, 6]
    # Verify that magnitude has changed times scalar value, angle hasn't changed
    initial_vector_magnitude = sqrt(initial.components[0]**2 + initial.components[1]**2 +
        initial.components[2]**2)
    result_vector_magnitude = sqrt(result[0]**2 + result[1]**2 + result[2]**2)
    assert result_vector_magnitude / initial_vector_magnitude == 2
    initial_vector_angle = atan2(initial.components[1], initial.components[0])
    result_vector_angle = atan2(result[1], result[0])
    assert initial_vector_angle == result_vector_angle

    assert scale_vector(1, Vector([2, 1])).components == [2, 1]
    assert scale_vector(0.1, Vector([1, 4])).components == [0.1, 0.4]
    assert scale_vector(2, Vector([1])).components == [2]
    assert len(scale_vector(2, Vector([])).components) == 0
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


def test_cylindrical_scale_vector(test_args: Args) -> None:
    vector1 = Vector([1, 2], test_args.C)
    Cy = coordinates_transform(test_args.C, CoordinateSystem.System.CYLINDRICAL)
    vector1_polar = vector1.rebase(Cy)
    scaled_vector1_polar = scale_vector(5, vector1_polar)
    scaled_vector1 = scaled_vector1_polar.rebase(test_args.C)

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
    vector2_polar = vector2.rebase(Cy)
    scaled_vector2_polar = scale_vector(5, vector2_polar)
    scaled_vector2 = scaled_vector2_polar.rebase(test_args.C)

    # Verify that magnitude has changed times scalar value, angle hasn't changed
    initial_vector2_magnitude = sqrt(vector2.components[0]**2 + vector2.components[1]**2 +
        vector2.components[2]**2)
    result_vector2_magnitude = sqrt(scaled_vector2.components[0]**2 +
        scaled_vector2.components[1]**2 + scaled_vector2.components[2]**2)
    assert result_vector2_magnitude / initial_vector2_magnitude == 5
    # It's hard to verify angle so only validate magnitude

    assert scale_vector(5, Vector([1, 2], Cy)).components == [5, 2]
    assert scale_vector(5, Vector([1, 2, 3], Cy)).components == [5, 2, 15]


def test_spherical_scale_vector(test_args: Args) -> None:
    vector1 = Vector([1, 2], test_args.C)
    Cs = coordinates_transform(test_args.C, CoordinateSystem.System.SPHERICAL)
    vector1_polar = vector1.rebase(Cs)
    scaled_vector1_polar = scale_vector(5, vector1_polar)
    scaled_vector1 = scaled_vector1_polar.rebase(test_args.C)

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
    vector2_polar = vector2.rebase(Cs)
    scaled_vector2_polar = scale_vector(5, vector2_polar)
    scaled_vector2 = scaled_vector2_polar.rebase(test_args.C)

    # Verify that magnitude has changed times scalar value, angle hasn't changed
    initial_vector2_magnitude = sqrt(vector2.components[0]**2 + vector2.components[1]**2 +
        vector2.components[2]**2)
    result_vector2_magnitude = sqrt(scaled_vector2.components[0]**2 +
        scaled_vector2.components[1]**2 + scaled_vector2.components[2]**2)
    assert result_vector2_magnitude / initial_vector2_magnitude == 5
    # It's hard to verify angle so only validate magnitude

    assert scale_vector(5, Vector([1, 2], Cs)).components == [5, 2]
    assert scale_vector(5, Vector([1, 2, 3], Cs)).components == [5, 2, 3]


def test_basic_dot_product(test_args: Args) -> None:
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
    initial = Vector([2, 3])
    initial_vector_angle = atan2(initial.components[1], initial.components[0])
    result_vector_angle = initial_vector_angle + pi / 2
    orthogonal_unit_vector = Vector([cos(result_vector_angle), sin(result_vector_angle)])
    assert dot_vectors(initial, orthogonal_unit_vector) == 0
    # Verify that parallel vectors have dot product equal to their magnitudes multiplication
    first_vector = Vector([1, 2])
    second_vector = Vector([2, 4])
    first_vector_magnitude = sqrt(first_vector.components[0]**2 + first_vector.components[1]**2)
    second_vector_magnitude = sqrt(second_vector.components[0]**2 + second_vector.components[1]**2)
    assert dot_vectors(first_vector,
        second_vector) == first_vector_magnitude * second_vector_magnitude


def test_invalid_dot_product(test_args: Args) -> None:
    with raises(ValueError):
        dot_vectors(Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
            Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]))
    with raises(ValueError):
        dot_vectors(Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
            Vector([test_args.C.coord_system.x * test_args.C.coord_system.j, 2]))
    C1 = CoordinateSystem()
    with raises(TypeError):
        dot_vectors(Vector([1, 2], test_args.C), Vector([1, 2], C1))


def test_non_cartesian_dot_product() -> None:
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    assert dot_vectors(Vector([1, 0], C1), Vector([1, pi / 2], C1)) == 0
    assert dot_vectors(Vector([1, 0], C1), Vector([1, pi / 4], C1)) == (sqrt(2) / 2)
    assert dot_vectors(Vector([1, 0, 0], C1), Vector([1, pi / 2, 0], C1)) == 0
    assert dot_vectors(Vector([1, 0, 0], C1), Vector([0, 0, 1], C1)) == 0

    # When r = 0 and z > 0, vector is directed up and is immune to rotations,
    # therefore rotated vector is parallel to original
    assert dot_vectors(Vector([0, 0, 1], C1), Vector([0, pi / 2, 1], C1)) == 1

    C2 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)
    assert dot_vectors(Vector([1, 0], C2), Vector([1, 0, pi / 2], C2)) == 0
    assert dot_vectors(Vector([1, 0], C2), Vector([1, 0, pi / 4], C2)) == (sqrt(2) / 2)
    assert dot_vectors(Vector([1, 0, 0], C2), Vector([1, 0, pi / 2], C2)) == 0
    assert dot_vectors(Vector([1, 0, pi / 2], C2), Vector([1, pi / 2, pi / 2], C2)) == 0

    # When phi angle = 0, vector is directed up and is immune to rotations,
    # therefore rotated vector is parallel to original
    assert dot_vectors(Vector([1, 0, 0], C2), Vector([1, pi / 2, 0], C2)) == 1


def test_basic_magnitude_vectors(test_args: Args) -> None:
    assert vector_magnitude(Vector([1, 2])) == sqrt(5)
    assert vector_magnitude(Vector([])) == 0
    assert vector_magnitude(Vector([1, 2, 3])) == sqrt(14)
    assert vector_magnitude(Vector([test_args.C.coord_system.x,
        2])) == sqrt(test_args.C.coord_system.x**2 + 4)


def test_non_cartesian_magnitude_vectors() -> None:
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    assert vector_magnitude(Vector([3, 0], C1)) == 3
    assert vector_magnitude(Vector([3, 0, 2], C1)) == sqrt(13)
    assert vector_magnitude(Vector([3, pi / 4, 2], C1)) == sqrt(13)

    C2 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)
    assert vector_magnitude(Vector([3, 0], C2)) == 3
    assert vector_magnitude(Vector([3, 2, 0], C2)) == 3
    assert vector_magnitude(Vector([3, 2, pi / 4], C2)) == 3


def test_basic_cross_product(test_args: Args) -> None:
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


def test_invalid_cross_product(test_args: Args) -> None:
    with raises(TypeError):
        cross_cartesian_vectors(
            Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
            Vector([2, test_args.C.coord_system.x * test_args.C.coord_system.i]))
    with raises(TypeError):
        cross_cartesian_vectors(
            Vector([test_args.C.coord_system.x * test_args.C.coord_system.i, 2]),
            Vector([2, test_args.C.coord_system.x * test_args.C.coord_system.j]))
    with raises(ValueError):
        cross_cartesian_vectors(Vector([1, 2, 3, 4]), Vector([1, 2]))
    C1 = CoordinateSystem()
    with raises(TypeError):
        cross_cartesian_vectors(Vector([1, 2], test_args.C), Vector([1, 2], C1))
    with raises(ValueError):
        cross_cartesian_vectors(Vector([1, 2, 3, 4], test_args.C), Vector([1, 2], test_args.C))
    # non-cartesian cross product is not supported
    C2 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    with raises(ValueError):
        cross_cartesian_vectors(Vector([1, 2], C2), Vector([1, 2], C2))
    C3 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)
    with raises(ValueError):
        cross_cartesian_vectors(Vector([1, 2], C3), Vector([1, 2], C3))


def test_rebased_cross_product(test_args: Args) -> None:
    Cy = coordinates_transform(test_args.C, CoordinateSystem.System.CYLINDRICAL)
    # First vector has r = 5 and theta = 0
    vector1 = Vector([5, 0, 1], Cy)
    # Second vector has r = 5 and theta = pi/2
    vector2 = Vector([5, pi / 2, 1], Cy)
    # Their product should be a vector with r = sqrt(50), theta = -3 * pi/4, z = 25
    vector1_cartesian = vector1.rebase(test_args.C)
    vector2_cartesian = vector2.rebase(test_args.C)
    vector_cartesian_cross = cross_cartesian_vectors(vector1_cartesian, vector2_cartesian)
    vector_cross = vector_cartesian_cross.rebase(Cy)
    assert vector_cross.components == [sqrt(50), -3 * pi / 4, 25]


def test_basic_unit_vector() -> None:
    assert vector_unit(Vector([1, 2])).components == [1 / sqrt(5), 2 / sqrt(5)]
    assert vector_unit(Vector([])).components == []
    assert vector_unit(Vector([1, 2, 3])).components == [1 / sqrt(14), 2 / sqrt(14), 3 / sqrt(14)]
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    assert vector_unit(Vector([3, 0], C1)).components == [1, 0]
    assert vector_unit(Vector([3, 0, 2], C1)).components == [3 / sqrt(13), 0, 2 / sqrt(13)]
    assert vector_unit(Vector([3, pi / 4, 2],
        C1)).components == [3 / sqrt(13), pi / 4, 2 / sqrt(13)]
    C2 = CoordinateSystem(CoordinateSystem.System.SPHERICAL)
    assert vector_unit(Vector([3, 0], C2)).components == [1, 0]
    assert vector_unit(Vector([3, 2, 0], C2)).components == [1, 2, 0]
    assert vector_unit(Vector([3, 2, pi / 4], C2)).components == [1, 2, pi / 4]


def test_basic_equal_quantity_vectors() -> None:
    Q1 = Quantity(1)
    Q2 = Quantity(2)
    assert equal_vectors(QuantityVector([Q1, Q2]), QuantityVector([Q1, Q2]))
    assert not equal_vectors(QuantityVector([Q1, Q2]), QuantityVector([Q2, Q1]))
    L1 = Quantity(1 * units.meter)
    L2 = Quantity(2 * units.meter)
    assert equal_vectors(QuantityVector([L1, L2]), QuantityVector([L1, L2]))
    assert not equal_vectors(QuantityVector([L1, L2]), QuantityVector([Q1, Q2]))


def test_basic_add_quantity_vectors() -> None:
    Q1 = Quantity(1)
    Q2 = Quantity(2)
    result = add_cartesian_quantity_vectors(QuantityVector([Q1, Q2]), QuantityVector([Q1, Q2]))
    assert [result.components[0].scale_factor, result.components[1].scale_factor] == [2, 4]
    result = add_cartesian_quantity_vectors(QuantityVector([Q1, Q2]), QuantityVector([Q2, Q1]))
    assert [result.components[0].scale_factor, result.components[1].scale_factor] == [3, 3]
    L1 = Quantity(1 * units.meter)
    L2 = Quantity(2 * units.meter)
    result = add_cartesian_quantity_vectors(QuantityVector([L1, L2]), QuantityVector([L1, L2]))
    assert [result.components[0].scale_factor, result.components[1].scale_factor] == [2, 4]
    assert result.dimension == units.length
    result = add_cartesian_quantity_vectors(QuantityVector([L1, L2]), QuantityVector([L2, L1]))
    assert [result.components[0].scale_factor, result.components[1].scale_factor] == [3, 3]
    assert result.dimension == units.length


def test_basic_scale_quantity_vectors() -> None:
    Q1 = Quantity(1)
    Q2 = Quantity(2)
    initial = QuantityVector([Q1, Q2])
    result = scale_quantity_vector(Quantity(2), initial)
    assert [result.components[0].scale_factor, result.components[1].scale_factor] == [2, 4]
    assert result.dimension == dimensionless
    L1 = Quantity(1 * units.meter)
    L2 = Quantity(2 * units.meter)
    initial = QuantityVector([L1, L2])
    result = scale_quantity_vector(Quantity(2), initial)
    assert [result.components[0].scale_factor, result.components[1].scale_factor] == [2, 4]
    assert result.dimension == units.length
    result = scale_quantity_vector(L2, initial)
    assert [result.components[0].scale_factor, result.components[1].scale_factor] == [2, 4]
    assert result.dimension == units.length**2


def test_basic_dot_quantity_product() -> None:
    Q1 = Quantity(1)
    Q2 = Quantity(2)
    assert dot_quantity_vectors(QuantityVector([Q1, Q2]), QuantityVector([Q1,
        Q2])).scale_factor == 5
    L1 = Quantity(1 * units.meter)
    L2 = Quantity(2 * units.meter)
    distance1 = QuantityVector([L1, L2])
    distance2 = QuantityVector([L2, L1])
    area = dot_quantity_vectors(distance1, distance2)
    assert area.scale_factor == 4
    SI.get_dimension_system().equivalent_dims(area.dimension, units.area)


def test_basic_quantity_magnitude() -> None:
    Q1 = Quantity(1)
    Q2 = Quantity(2)
    magnitude = quantity_vector_magnitude(QuantityVector([Q1, Q2]))
    assert magnitude.scale_factor == sqrt(5)
    assert magnitude.dimension == dimensionless
    L1 = Quantity(1 * units.meter)
    L2 = Quantity(2 * units.meter)
    initial = QuantityVector([L1, L2])
    magnitude = quantity_vector_magnitude(initial)
    assert magnitude.scale_factor == sqrt(5)
    SI.get_dimension_system().equivalent_dims(magnitude.dimension, units.length)


def test_basic_quantity_cross_product() -> None:
    Q1 = Quantity(1)
    Q2 = Quantity(2)
    crossed = cross_cartesian_quantity_vectors(QuantityVector([Q1, Q2]), QuantityVector([Q1, Q2]))
    assert [
        crossed.components[0].scale_factor, crossed.components[1].scale_factor,
        crossed.components[2].scale_factor
    ] == [0, 0, 0]
    crossed = cross_cartesian_quantity_vectors(QuantityVector([Q1, Q2]), QuantityVector([Q2, Q1]))
    assert [
        crossed.components[0].scale_factor, crossed.components[1].scale_factor,
        crossed.components[2].scale_factor
    ] == [0, 0, -3]
    L1 = Quantity(1 * units.meter)
    L2 = Quantity(2 * units.meter)
    distance1 = QuantityVector([L1, L2])
    distance2 = QuantityVector([L2, L1])
    cross_area_vector = cross_cartesian_quantity_vectors(distance1, distance2)
    assert [
        cross_area_vector.components[0].scale_factor, cross_area_vector.components[1].scale_factor,
        cross_area_vector.components[2].scale_factor
    ] == [0, 0, -3]
    SI.get_dimension_system().equivalent_dims(cross_area_vector.dimension, units.area)


def test_basic_quantity_unit_vector() -> None:
    Q1 = Quantity(1)
    Q2 = Quantity(2)
    unit_vector = quantity_vector_unit(QuantityVector([Q1, Q2]))
    assert [unit_vector.components[0].scale_factor,
        unit_vector.components[1].scale_factor] == [1 / sqrt(5), 2 / sqrt(5)]
    L1 = Quantity(1 * units.meter)
    L2 = Quantity(2 * units.meter)
    initial = QuantityVector([L1, L2])
    unit_vector = quantity_vector_unit(initial)
    assert [unit_vector.components[0].scale_factor,
        unit_vector.components[1].scale_factor] == [1 / sqrt(5), 2 / sqrt(5)]
    SI.get_dimension_system().equivalent_dims(unit_vector.components[0].dimension, units.length)
    SI.get_dimension_system().equivalent_dims(unit_vector.components[1].dimension, units.length)
