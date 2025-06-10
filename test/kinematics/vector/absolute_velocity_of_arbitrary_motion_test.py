from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.kinematics.vector import absolute_velocity_of_arbitrary_motion as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors
from symplyphysics.core.experimental.solvers import solve_for_vector

Args = namedtuple("Args", "vabs vrel vtr")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    unit = units.meter / units.second

    vabs = QuantityCoordinateVector([-1 * unit, 3 * unit, -4 * unit], CARTESIAN)
    vrel = QuantityCoordinateVector([0, 1 * unit, -1 * unit], CARTESIAN)
    vtr = QuantityCoordinateVector([-1 * unit, 2 * unit, -3 * unit], CARTESIAN)
    return Args(vabs=vabs, vrel=vrel, vtr=vtr)


def test_absolute_law(test_args: Args) -> None:
    result = law.calculate_absolute_velocity(test_args.vrel, test_args.vtr)
    assert_equal_vectors(result, test_args.vabs)


def test_relative_law(test_args: Args) -> None:
    expr = solve_for_vector(law.law, law.relative_velocity)

    result = expr.subs({
        law.absolute_velocity: test_args.vabs,
        law.transfer_velocity: test_args.vtr,
    })

    result = QuantityCoordinateVector.from_expr(result)

    assert_equal_vectors(result, test_args.vrel)


def test_transfer_law(test_args: Args) -> None:
    expr = solve_for_vector(law.law, law.transfer_velocity)

    result = expr.subs({
        law.absolute_velocity: test_args.vabs,
        law.relative_velocity: test_args.vrel,
    })

    result = QuantityCoordinateVector.from_expr(result)

    assert_equal_vectors(result, test_args.vtr)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityCoordinateVector([Quantity(1 * units.coulomb), 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_absolute_velocity(vb_vector, test_args.vtr)
    with raises(errors.UnitsError):
        law.calculate_absolute_velocity(test_args.vrel, vb_vector)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(ValueError):
        law.calculate_absolute_velocity(vb_scalar, test_args.vtr)
    with raises(ValueError):
        law.calculate_absolute_velocity(test_args.vrel, vb_scalar)

    with raises(TypeError):
        law.calculate_absolute_velocity(100, test_args.vtr)
    with raises(TypeError):
        law.calculate_absolute_velocity([100], test_args.vtr)
    with raises(TypeError):
        law.calculate_absolute_velocity(test_args.vrel, 100)
    with raises(TypeError):
        law.calculate_absolute_velocity(test_args.vrel, [100])
