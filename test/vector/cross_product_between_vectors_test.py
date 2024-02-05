from collections import namedtuple
from pytest import approx, fixture, raises, mark
from symplyphysics import Vector
from symplyphysics.laws.vector import (
    cross_product_between_vectors as cross_law,
)

# Description
## Given two vectors, a = (1, 1, 0) and b = (0, 2, -1), the norm of the cross product
## between them is approximately 2.45


@fixture(name="test_args")
def test_args_fixture():
    a = Vector([1, 1, 0])
    b = Vector([0, 2, -1])
    c = Vector([0, 0, 4])  # perpendicular to a
    d = Vector([2, 2, 0])  # parallel to a
    Args = namedtuple("Args", "a b c d")
    return Args(a=a, b=b, c=c, d=d)


@mark.skip("")
def test_basic_law(test_args):
    result = cross_law.cross_product_law(test_args.a, test_args.b)
    assert result.evalf(3) == approx(2.45, 1e-3)


@mark.skip("")
def test_basic_law_perpendicular_vectors(test_args):
    result = cross_law.cross_product_law(test_args.a, test_args.c)
    assert result.evalf(3) == approx(5.66, 1e-3)


@mark.skip("")
def test_basic_law_parallel_vectors(test_args):
    result = cross_law.cross_product_law(test_args.a, test_args.d)
    assert result.evalf(3) == approx(0, abs=1e-14)


@mark.skip("")
def test_bad_dimensions_count(test_args):
    a_bad = Vector([1, 2, 3, 4, 5])
    b_bad = Vector([-1, -2, -3, -4, 5])
    with raises(ValueError):
        cross_law.cross_product_law(a_bad, b_bad)
