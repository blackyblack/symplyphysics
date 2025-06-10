from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.gravity.vector import acceleration_due_to_gravity_via_gravity_force_and_centripetal_acceleration as law

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector, CARTESIAN
from symplyphysics.core.experimental.approx import assert_equal_vectors
from symplyphysics.core.experimental.solvers import solve_for_vector

Args = namedtuple("Args", "f a m g")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f = QuantityCoordinateVector([
        -2.43e2 * units.newton,
        -4.21e2 * units.newton,
        -4.86e2 * units.newton,
    ], CARTESIAN)

    a_unit = units.meter / units.second**2
    a = QuantityCoordinateVector([
        -1.2e-2 * a_unit,
        -2.07e-2 * a_unit,
        0,
    ], CARTESIAN)

    m = Quantity(70 * units.kilogram)

    g = QuantityCoordinateVector([
        -3.46 * a_unit,
        -5.99 * a_unit,
        -6.94 * a_unit,
    ], CARTESIAN)

    return Args(f=f, a=a, m=m, g=g)


def test_law(test_args: Args) -> None:
    result = law.calculate_acceleraton_due_to_gravity(test_args.f, test_args.a, test_args.m)
    assert_equal_vectors(result, test_args.g)


def test_gravity_law(test_args: Args) -> None:
    result = solve_for_vector(law.law, law.gravity_force).subs({
        law.mass: test_args.m,
        law.acceleration_due_to_gravity: test_args.g,
        law.centripetal_acceleration: test_args.a,
    })
    result = QuantityCoordinateVector.from_expr(result)

    assert_equal_vectors(result, test_args.f)


def test_bad_force(test_args: Args) -> None:
    fb_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_acceleraton_due_to_gravity(fb_vector, test_args.a, test_args.m)

    fb_scalar = Quantity(units.newton)
    with raises(ValueError):
        law.calculate_acceleraton_due_to_gravity(fb_scalar, test_args.a, test_args.m)

    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity(100, test_args.a, test_args.m)
    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity([100], test_args.a, test_args.m)


def test_bad_acceleration(test_args: Args) -> None:
    ab_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, ab_vector, test_args.m)

    ab_scalar = Quantity(units.meter / units.second**2)
    with raises(ValueError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, ab_scalar, test_args.m)

    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, 100, test_args.m)
    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, [100], test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, test_args.a, mb)
    with raises(TypeError):
        law.calculate_acceleraton_due_to_gravity(test_args.f, test_args.a, 100)
