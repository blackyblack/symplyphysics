from collections import namedtuple
from pytest import fixture
from symplyphysics import assert_equal
from symplyphysics.core.coordinate_systems import CARTESIAN, CoordinateVector
from symplyphysics.core.vectors import VectorCross, VectorNorm
from symplyphysics.laws.mathematics.geometry import cross_product_is_proportional_to_sine_between_vectors as sine_law

# Description
## A force F = (1, 2, -1) N is acting on a point at position r = (0, 1, -2) m. Therefore
## the torque of the force acting on the point is tau = cross(r, F) = (3, -2, -1) N*m.
## The sine of the angle between F and r is 0.638.

Args = namedtuple("Args", "r_norm f_norm tau_norm")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f = CoordinateVector([1.0, 2.0, -1.0], CARTESIAN)
    r = CoordinateVector([0.0, 1.0, -2.0], CARTESIAN)
    tau = VectorCross(r, f)
    r_norm = VectorNorm(r)
    f_norm = VectorNorm(f)
    tau_norm = VectorNorm(tau)
    return Args(r_norm=r_norm, f_norm=f_norm, tau_norm=tau_norm)


def test_law(test_args: Args) -> None:
    result = sine_law.calculate_sine_between_vectors(test_args.tau_norm, test_args.r_norm,
        test_args.f_norm)
    assert_equal(result, 0.683)


def test_law_perpendicular(test_args: Args) -> None:
    tau_perpendicular_norm = test_args.r_norm * test_args.f_norm
    result = sine_law.calculate_sine_between_vectors(tau_perpendicular_norm, test_args.r_norm,
        test_args.f_norm)
    assert_equal(result, 1.0)


def test_law_parallel(test_args: Args) -> None:
    tau_parallel_norm = 0.0
    result = sine_law.calculate_sine_between_vectors(tau_parallel_norm, test_args.r_norm,
        test_args.f_norm)
    assert_equal(result, 0.0)
