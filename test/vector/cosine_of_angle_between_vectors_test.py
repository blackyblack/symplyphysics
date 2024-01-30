from collections import namedtuple
from pytest import approx, fixture
from sympy import pi
from symplyphysics import Vector, vector_magnitude, dot_vectors
from symplyphysics.laws.vector import cosine_of_angle_between_vectors as cosine_law

# Description
## Given two vectors a = (-1, 2, 0.5) and b = (0, 3, 0), the cosine of the angle
## between the two is close to 0.873, which corresponds to an angle of 29.2 degrees.


@fixture(name="test_args")
def test_args_fixture():
    a = Vector([-1, 2, 0.5])
    b = Vector([0, 3, 0])

    a_dot_b = dot_vectors(a, b)
    a_norm = vector_magnitude(a)
    b_norm = vector_magnitude(b)

    Args = namedtuple("Args", "a b a_dot_b a_norm b_norm")
    return Args(a=a, b=b, a_dot_b=a_dot_b, a_norm=a_norm, b_norm=b_norm)


def test_vector_law(test_args):
    result = cosine_law.cosine_law(test_args.a, test_args.b)
    assert result.evalf(3) == approx(0.873, 1e-3)


def test_basic_law(test_args):
    result_rad = cosine_law.calculate_angle_between_vectors(
        test_args.a_dot_b, test_args.a_norm, test_args.b_norm
    )
    result_deg = (result_rad * 180 / pi).evalf(3)
    assert result_deg == approx(29.2, 1e-3)


def test_basic_law_perpendicular(test_args):
    result = cosine_law.calculate_angle_between_vectors(0, test_args.a_norm, test_args.b_norm)
    assert result == approx(pi / 2, 1e-3)


def test_basic_law_parallel(test_args):
    result = cosine_law.calculate_angle_between_vectors(
        test_args.a_norm * test_args.b_norm,
        test_args.a_norm,
        test_args.b_norm,
    )
    assert result == approx(0, abs=1e-14)
