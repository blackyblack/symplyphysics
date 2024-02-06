from collections import namedtuple
from pytest import approx, fixture
from symplyphysics import (
    Quantity,
    QuantityVector,
    dimensionless,
    SI,
    units,
    convert_to,
)
from symplyphysics.laws.geometry.vector import (
    dot_product_is_proportional_to_cosine_between_vectors as dot_product_law,
)

# Description
## A particle with position vector (0.0, 0.3, 2.3) m is moving under a force (1.0, 2.0, -1.5) N.
## The cosine between the force vector and the position vector amounts to about -0.456.


@fixture(name="test_args")
def test_args_fixture():
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
    Args = namedtuple("Args", "f r")
    return Args(f=f, r=r)


def test_law(test_args):
    result = dot_product_law.calculate_cosine_between_vectors(test_args.f, test_args.r)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, dimensionless)
    result_value = convert_to(result, dimensionless).evalf(3)
    assert result_value == approx(-0.456, 1e-3)


def test_law_perpendicular(test_args):
    f_perpendicular = QuantityVector([
        Quantity(1.0 * units.newton),
        Quantity(0.0 * units.newton),
        Quantity(0.0 * units.newton),
    ])
    result = dot_product_law.calculate_cosine_between_vectors(f_perpendicular, test_args.r)
    result_value = convert_to(result, dimensionless).evalf(3)
    assert result_value == approx(0.0, abs=1e-15)


def test_law_parallel(test_args):
    f_parallel = QuantityVector([
        Quantity(0.0 * units.newton),
        Quantity(0.3 * units.newton),
        Quantity(2.3 * units.newton),
    ])
    result = dot_product_law.calculate_cosine_between_vectors(f_parallel, test_args.r)
    result_value = convert_to(result, dimensionless).evalf(3)
    assert result_value == approx(1.0, 1e-3)
