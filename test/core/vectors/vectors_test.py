from collections import namedtuple
from pytest import fixture, raises
from sympy import atan, pi, sqrt, symbols, sin, cos
from sympy.vector import Vector as SympyVector, express
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem, coordinates_rotate, coordinates_transform
from symplyphysics.core.vectors.vectors import Vector, sympy_vector_from_vector, vector_from_sympy_vector, vector_rebase


@fixture(name="test_args")
def test_args_fixture():
    C = CoordinateSystem()
    Args = namedtuple("Args", ["C"])
    return Args(C=C)


# Test Vector constructor


def test_basic_vector():
    vector = Vector([1, 2])
    assert vector.components == [1, 2]
    assert vector.coordinate_system is None


def test_coord_sys_vector(test_args):
    vector = Vector([1, 2], test_args.C)
    assert vector.components == [1, 2]
    assert vector.coordinate_system == test_args.C


def test_empty_vector():
    vector = Vector([])
    assert len(vector.components) == 0
    assert vector.coordinate_system is None


# Test vector_from_sympy_vector()


def test_basic_sympy_to_array_conversion(test_args):
    vector = vector_from_sympy_vector(test_args.C.coord_system.i + 2 * test_args.C.coord_system.j,
        test_args.C)
    assert vector.components == [1, 2, 0]
    assert vector.coordinate_system == test_args.C


def test_order_sympy_to_array_conversion(test_args):
    vector = vector_from_sympy_vector(2 * test_args.C.coord_system.j + test_args.C.coord_system.i,
        test_args.C)
    assert vector.components == [1, 2, 0]


def test_skip_dimension_sympy_to_array_conversion(test_args):
    vector = vector_from_sympy_vector(test_args.C.coord_system.i + 2 * test_args.C.coord_system.k,
        test_args.C)
    assert vector.components == [1, 0, 2]


def test_empty_sympy_to_array_conversion(test_args):
    vector = vector_from_sympy_vector(SympyVector.zero, test_args.C)
    assert len(vector.components) == 0
    assert vector.coordinate_system == test_args.C


# Does not support non SymPy Vectors
def test_only_scalar_sympy_to_array_conversion(test_args):
    with raises(AttributeError):
        vector_from_sympy_vector(test_args.C.coord_system.x, test_args.C)
    x1 = symbols("x1")
    with raises(AttributeError):
        vector_from_sympy_vector(x1, test_args.C)


def test_free_variable_sympy_to_array_conversion(test_args):
    x1 = symbols("x1")
    vector = vector_from_sympy_vector(test_args.C.coord_system.i * x1, test_args.C)
    assert vector.components == [x1, 0, 0]
    assert vector.coordinate_system == test_args.C


def test_non_cartesian_array_to_sympy_conversion():
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    vector = vector_from_sympy_vector(C1.coord_system.i + 2 * C1.coord_system.j, C1)
    assert vector.components == [1, 2, 0]
    assert vector.coordinate_system == C1


def test_rotate_coordinates_array_to_sympy_conversion(test_args):
    sympy_vector = test_args.C.coord_system.i + test_args.C.coord_system.j
    vector = vector_from_sympy_vector(sympy_vector, test_args.C)
    assert vector.components == [1, 1, 0]
    theta = symbols("theta")
    B = coordinates_rotate(test_args.C, theta, test_args.C.coord_system.k)
    sympy_transformed_vector = express(sympy_vector, B.coord_system)
    assert sympy_transformed_vector == ((sin(theta) + cos(theta)) * B.coord_system.i +
        (-sin(theta) + cos(theta)) * B.coord_system.j)
    transformed_vector = vector_from_sympy_vector(sympy_transformed_vector, B)
    assert transformed_vector.components == [sin(theta) + cos(theta), -sin(theta) + cos(theta), 0]
    assert transformed_vector.coordinate_system == B
    # Do the same via vector_rebase
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.components == [sin(theta) + cos(theta), -sin(theta) + cos(theta), 0]


# Test sympy_vector_from_vector()


def test_multiple_coord_systems_sympy_to_array_conversion(test_args):
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    with raises(TypeError):
        vector_from_sympy_vector(test_args.C.coord_system.i + 2 * C1.coord_system.k, test_args.C)
    with raises(TypeError):
        vector_from_sympy_vector(
            test_args.C.coord_system.i + C1.coord_system.r * test_args.C.coord_system.k, C1)


def test_basic_array_to_sympy_conversion(test_args):
    sympy_vector = sympy_vector_from_vector(Vector([1, 2], test_args.C))
    assert sympy_vector == test_args.C.coord_system.i + 2 * test_args.C.coord_system.j


def test_skip_dimension_array_to_sympy_conversion(test_args):
    sympy_vector = sympy_vector_from_vector(Vector([1, 0, 2], test_args.C))
    assert sympy_vector == test_args.C.coord_system.i + 2 * test_args.C.coord_system.k


def test_4d_array_to_sympy_conversion(test_args):
    sympy_vector = sympy_vector_from_vector(Vector([1, 0, 2, 5], test_args.C))
    assert sympy_vector == test_args.C.coord_system.i + 2 * test_args.C.coord_system.k


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
    B = coordinates_rotate(test_args.C, theta, test_args.C.coord_system.k)
    sympy_vector = sympy_vector_from_vector(Vector([1, 2], B))
    assert sympy_vector == B.coord_system.i + 2 * B.coord_system.j
    transformed_vector = express(sympy_vector, test_args.C.coord_system)
    assert transformed_vector == ((-2 * sin(theta) + cos(theta)) * test_args.C.coord_system.i +
        (sin(theta) + 2 * cos(theta)) * test_args.C.coord_system.j)
    # Do the same via vector_rebase
    vector_rebased = vector_rebase(Vector([1, 2], B), test_args.C)
    assert vector_rebased.components == [
        -2 * sin(theta) + cos(theta), sin(theta) + 2 * cos(theta), 0
    ]


# Test vector_rebase()


def test_basic_vector_rebase(test_args):
    vector = Vector([test_args.C.coord_system.x, test_args.C.coord_system.y], test_args.C)

    # B is located at [1, 2] origin instead of [0, 0] of test_args.C
    Bi = test_args.C.coord_system.locate_new(
        "B", test_args.C.coord_system.i + 2 * test_args.C.coord_system.j)
    B = CoordinateSystem(test_args.C.coord_system_type, Bi)
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coordinate_system == B
    # Original field is not changed
    assert vector.coordinate_system == test_args.C
    assert vector_rebased.components == [B.coord_system.x + 1, B.coord_system.y + 2, 0]


# Simple numbers in vector are not scalars - they cannot be properly
# rebased to new coordinate system.
def test_plain_vector_rebase(test_args):
    vector = Vector([1, 2], test_args.C)

    # B is located at [1, 2] origin instead of [0, 0] of test_args.C
    Bi = test_args.C.coord_system.locate_new(
        "B", test_args.C.coord_system.i + 2 * test_args.C.coord_system.j)
    B = CoordinateSystem(test_args.C.coord_system_type, Bi)
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coordinate_system == B
    # Original field is not changed
    assert vector.coordinate_system == test_args.C
    assert vector_rebased.components == [1, 2, 0]


# Simple parameters in vector are not scalars - they cannot be properly
# rebased to new coordinate system.
def test_parameters_vector_rebase(test_args):
    parameter = symbols("parameter")
    vector = Vector([parameter, parameter], test_args.C)

    # B is located at [1, 2] origin instead of [0, 0] of test_args.C
    Bi = test_args.C.coord_system.locate_new(
        "B", test_args.C.coord_system.i + 2 * test_args.C.coord_system.j)
    B = CoordinateSystem(test_args.C.coord_system_type, Bi)
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coordinate_system == B
    # Original field is not changed
    assert vector.coordinate_system == test_args.C
    assert vector_rebased.components == [parameter, parameter, 0]


# Vector is not rebased if no original coordinate system was set.
def test_no_coord_system_vector_rebase(test_args):
    vector = Vector([test_args.C.coord_system.x, test_args.C.coord_system.y])
    assert vector.coordinate_system is None
    Bi = test_args.C.coord_system.locate_new(
        "B", test_args.C.coord_system.i + 2 * test_args.C.coord_system.j)
    B = CoordinateSystem(test_args.C.coord_system_type, Bi)
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coordinate_system == B
    assert vector_rebased.components == [test_args.C.coord_system.x, test_args.C.coord_system.y]


# Vector is not rebased if no target coordinate system was set.
def test_no_target_coord_system_vector_rebase(test_args):
    vector = Vector([test_args.C.coord_system.x, test_args.C.coord_system.y], test_args.C)
    assert vector.coordinate_system == test_args.C
    vector_rebased = vector_rebase(vector, None)
    assert vector_rebased.coordinate_system is None
    assert vector_rebased.components == [test_args.C.coord_system.x, test_args.C.coord_system.y]


# Rotation does not require vector defined with base scalars.
def test_rotate_vector_rebase(test_args):
    vector = Vector([1, 2], test_args.C)
    point = [1, 2]
    p1 = test_args.C.coord_system.origin.locate_new(
        "p1", point[0] * test_args.C.coord_system.i + point[1] * test_args.C.coord_system.j)
    p1_coordinates = p1.express_coordinates(test_args.C.coord_system)
    assert p1_coordinates[0] == point[0]
    assert p1_coordinates[1] == point[1]

    B = coordinates_rotate(test_args.C, pi / 4, test_args.C.coord_system.k)
    p1_coordinates_in_b = p1.express_coordinates(B.coord_system)
    assert p1_coordinates_in_b[0] != point[0]

    transformed_point = [p1_coordinates_in_b[0], p1_coordinates_in_b[1], 0]
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coordinate_system == B
    assert vector_rebased.components == [3 * sqrt(2) / 2, sqrt(2) / 2, 0]
    assert vector_rebased.components == transformed_point


# Test non-cartesian coordinate systems


def test_spherical_vector_create(test_args):
    vector = Vector([1, 2], test_args.C)
    B = coordinates_transform(test_args.C, CoordinateSystem.System.SPHERICAL)
    # vector should have r = sqrt(5) in polar coordinates
    # vector is in XY-plane, so theta angle should be pi/2
    # phi angle is atan(2/1)
    vector_rebased = vector_rebase(vector, B)
    assert vector_rebased.coordinate_system == B
    assert vector_rebased.components == [sqrt(5), pi / 2, atan(2)]
