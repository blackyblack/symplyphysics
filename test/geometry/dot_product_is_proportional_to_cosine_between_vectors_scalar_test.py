from collections import namedtuple
from pytest import fixture
from symplyphysics import (
    assert_equal,
    Quantity,
    QuantityVector,
    vector_magnitude,
    dot_vectors,
    units,
)
from symplyphysics.laws.geometry import dot_product_is_proportional_to_cosine_between_vectors as cosine_law

# Description
## A point with position vector of (0.0, 0.3, 2.3) m is moving under a force (1.0, 2.0, -1.5) N.
## The cosine between the two vectors amounts to about -0.456.

Args = namedtuple("Args", "f_norm r_norm f_dot_r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f = QuantityVector([
        Quantity(1.0 * units.newton),
        Quantity(2.0 * units.newton),
        Quantity(-1.5 * units.newton),
    ])
    r = QuantityVector([
        Quantity(0.0 * units.meter),
        Quantity(0.3 * units.meter),
        Quantity(2.3 * units.meter),
    ])

    f_norm = vector_magnitude(f)
    r_norm = vector_magnitude(r)
    f_dot_r = dot_vectors(f, r)
    return Args(f_norm, r_norm, f_dot_r)


def test_law(test_args: Args) -> None:
    result = cosine_law.calculate_cosine_between_vectors(
        test_args.f_dot_r,
        test_args.f_norm,
        test_args.r_norm,
    )
    assert_equal(result, -0.456)


def test_law_perpendicular(test_args: Args) -> None:
    f_perp_dot_r = 0.0
    result = cosine_law.calculate_cosine_between_vectors(
        f_perp_dot_r,
        test_args.f_norm,
        test_args.r_norm,
    )
    assert_equal(result, 0.0)


def test_law_parallel(test_args: Args) -> None:
    f_parallel_dot_r = test_args.f_norm * test_args.r_norm
    result = cosine_law.calculate_cosine_between_vectors(
        f_parallel_dot_r,
        test_args.f_norm,
        test_args.r_norm,
    )
    assert_equal(result, 1.0)
