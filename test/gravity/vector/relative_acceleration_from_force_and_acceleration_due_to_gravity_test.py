from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.gravity.vector import (
    relative_acceleration_from_force_and_acceleration_due_to_gravity as law)

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors
from symplyphysics.core.experimental.solvers import solve_for_vector

Args = namedtuple("Args", "g ac f m ar")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a_unit = units.meter / units.second**2
    g = QuantityCoordinateVector([0, 0, -9.81 * a_unit], CARTESIAN)
    ac = QuantityCoordinateVector([0, -0.1 * a_unit, 0.25 * a_unit], CARTESIAN)
    f = QuantityCoordinateVector([100 * units.newton, 100 * units.newton, 0], CARTESIAN)
    m = Quantity(100 * units.kilogram)
    ar = QuantityCoordinateVector([1.0 * a_unit, 1.1 * a_unit, -10.06 * a_unit], CARTESIAN)
    return Args(g=g, ac=ac, f=f, m=m, ar=ar)


def test_relative_acceleration_law(test_args: Args) -> None:
    result = law.calculate_acceleration(test_args.g, test_args.ac, test_args.f, test_args.m)
    assert_equal_vectors(result, test_args.ar)


def test_gravity_acceleration_law(test_args: Args) -> None:
    expr = solve_for_vector(law.law, law.acceleration_due_to_gravity)

    result = expr.subs({
        law.relative_acceleration: test_args.ar,
        law.coriolis_acceleration: test_args.ac,
        law.mass: test_args.m,
        law.force: test_args.f,
    })

    result = QuantityCoordinateVector.from_expr(result)

    assert_equal_vectors(result, test_args.g, absolute_tolerance=1e-10)


def test_coriolis_acceleration_law(test_args: Args) -> None:
    expr = solve_for_vector(law.law, law.coriolis_acceleration)

    result = expr.subs({
        law.relative_acceleration: test_args.ar,
        law.acceleration_due_to_gravity: test_args.g,
        law.mass: test_args.m,
        law.force: test_args.f,
    })

    result = QuantityCoordinateVector.from_expr(result)

    assert_equal_vectors(result, test_args.ac)


def test_force_law(test_args: Args) -> None:
    expr = solve_for_vector(law.law, law.force)

    result = expr.subs({
        law.relative_acceleration: test_args.ar,
        law.acceleration_due_to_gravity: test_args.g,
        law.coriolis_acceleration: test_args.ac,
        law.mass: test_args.m,
    })

    result = QuantityCoordinateVector.from_expr(result)

    assert_equal_vectors(result, test_args.f)


def test_bad_acceleration(test_args: Args) -> None:
    ab_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_acceleration(ab_vector, test_args.ac, test_args.f, test_args.m)
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.g, ab_vector, test_args.f, test_args.m)

    ab_scalar = units.meter / units.second**2
    with raises(ValueError):
        law.calculate_acceleration(ab_scalar, test_args.ac, test_args.f, test_args.m)
    with raises(ValueError):
        law.calculate_acceleration(test_args.g, ab_scalar, test_args.f, test_args.m)

    with raises(TypeError):
        law.calculate_acceleration(100, test_args.ac, test_args.f, test_args.m)
    with raises(TypeError):
        law.calculate_acceleration([100], test_args.ac, test_args.f, test_args.m)
    with raises(TypeError):
        law.calculate_acceleration(test_args.g, 100, test_args.f, test_args.m)
    with raises(TypeError):
        law.calculate_acceleration(test_args.g, [100], test_args.f, test_args.m)


def test_bad_force(test_args: Args) -> None:
    fb_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.g, test_args.ac, fb_vector, test_args.m)

    fb_scalar = units.newton
    with raises(ValueError):
        law.calculate_acceleration(test_args.g, test_args.ac, fb_scalar, test_args.m)

    with raises(TypeError):
        law.calculate_acceleration(test_args.g, test_args.ac, 100, test_args.m)
    with raises(TypeError):
        law.calculate_acceleration(test_args.g, test_args.ac, [100], test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.g, test_args.ac, test_args.f, mb)
    with raises(TypeError):
        law.calculate_acceleration(test_args.g, test_args.ac, test_args.f, 100)
