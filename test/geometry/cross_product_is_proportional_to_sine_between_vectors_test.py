from collections import namedtuple
from pytest import fixture
from symplyphysics import (
    assert_equal,
    Quantity,
    QuantityVector,
    vector_magnitude,
    cross_cartesian_vectors,
    units,
)
from symplyphysics.laws.geometry import (
    cross_product_is_proportional_to_sine_between_vectors as sine_law,)

# Description
## A force F = (1, 2, -1) N is acting on a point at position r = (0, 1, -2) m. Therefore
## the torque of the force acting on the point is tau = cross(r, F) = (3, -2, -1) N*m.
## The sine of the angle between F and r is 0.638.


@fixture(name="test_args")
def test_args_fixture():
    r = QuantityVector([
        Quantity(0.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(-2.0 * units.meter),
    ])
    F = QuantityVector([
        Quantity(1.0 * units.newton),
        Quantity(2.0 * units.newton),
        Quantity(-1.0 * units.newton),
    ])
    tau = cross_cartesian_vectors(r, F)

    r_norm = vector_magnitude(r)
    f_norm = vector_magnitude(F)
    tau_norm = vector_magnitude(tau)

    Args = namedtuple("Args", "r_norm f_norm tau_norm")
    return Args(r_norm=r_norm, f_norm=f_norm, tau_norm=tau_norm)


def test_law(test_args):
    result = sine_law.calculate_sine_between_vectors(test_args.tau_norm, test_args.r_norm,
        test_args.f_norm)
    assert_equal(result, 0.683)


def test_law_perpendicular(test_args):
    tau_perpendicular_norm = test_args.r_norm * test_args.f_norm
    result = sine_law.calculate_sine_between_vectors(tau_perpendicular_norm, test_args.r_norm,
        test_args.f_norm)
    assert_equal(result, 1.0)


def test_law_parallel(test_args):
    tau_parallel_norm = 0.0
    result = sine_law.calculate_sine_between_vectors(tau_parallel_norm, test_args.r_norm,
        test_args.f_norm)
    assert_equal(result, 0.0)
