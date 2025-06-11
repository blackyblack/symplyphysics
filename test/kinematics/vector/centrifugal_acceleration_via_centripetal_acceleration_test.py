from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units
from symplyphysics.laws.kinematics.vector import centrifugal_acceleration_via_centripetal_acceleration as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors
from symplyphysics.core.experimental.solvers import solve_for_vector

Args = namedtuple("Args", "acp acf")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a_unit = units.meter / units.second**2

    acp = QuantityCoordinateVector([2 * a_unit, -1 * a_unit, 3 * a_unit], CARTESIAN)
    acf = QuantityCoordinateVector([-2 * a_unit, 1 * a_unit, -3 * a_unit], CARTESIAN)

    return Args(acp=acp, acf=acf)


def test_centrifugal_law(test_args: Args) -> None:
    result = law.calculate_centrifugal_acceleration(test_args.acp)
    assert_equal_vectors(result, test_args.acf)


def test_centripetal_law(test_args: Args) -> None:
    result = solve_for_vector(law.law, law.centripetal_acceleration).subs(
        law.centrifugal_acceleration,
        test_args.acf,
    )

    result = QuantityCoordinateVector.from_expr(result)
    assert_equal_vectors(result, test_args.acp)


def test_bad_acceleration() -> None:
    ab_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_centrifugal_acceleration(ab_vector)

    ab_scalar = units.meter / units.second**2
    with raises(ValueError):
        law.calculate_centrifugal_acceleration(ab_scalar)

    with raises(TypeError):
        law.calculate_centrifugal_acceleration(100)
    with raises(TypeError):
        law.calculate_centrifugal_acceleration([100])
