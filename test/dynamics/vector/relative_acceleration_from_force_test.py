from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.dynamics.vector import relative_acceleration_from_force as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors
from symplyphysics.core.experimental.solvers import solve_for_vector

Args = namedtuple("Args", "m a_rel f a_cor a_tr")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(5 * units.kilogram)

    a_unit = units.meter / units.second**2
    a_rel = QuantityCoordinateVector([6 * a_unit, 1 * a_unit, -5 * a_unit], CARTESIAN)
    a_cor = QuantityCoordinateVector([0, 4 * a_unit, -2 * a_unit], CARTESIAN)
    a_tr = QuantityCoordinateVector([-2 * a_unit, 3 * a_unit, 0], CARTESIAN)

    f = QuantityCoordinateVector([20 * units.newton, 0, -15 * units.newton], CARTESIAN)

    return Args(m=m, a_rel=a_rel, f=f, a_cor=a_cor, a_tr=a_tr)


def test_relative_law(test_args: Args) -> None:
    result = law.calculate_relative_acceleration(test_args.m, test_args.f, test_args.a_cor,
        test_args.a_tr)
    assert_equal_vectors(result, test_args.a_rel)


def test_force_law(test_args: Args) -> None:
    force_expr = solve_for_vector(law.law, law.force)

    force_subs = force_expr.subs({
        law.relative_acceleration: test_args.a_rel,
        law.mass: test_args.m,
        law.coriolis_acceleration: test_args.a_cor,
        law.translation_acceleration: test_args.a_tr,
    })

    force_value = QuantityCoordinateVector.from_expr(force_subs)

    assert_equal_vectors(force_value, test_args.f)


def test_coriolis_law(test_args: Args) -> None:
    acceleration_expr = solve_for_vector(law.law, law.coriolis_acceleration)

    acceleration_subs = acceleration_expr.subs({
        law.relative_acceleration: test_args.a_rel,
        law.force: test_args.f,
        law.mass: test_args.m,
        law.translation_acceleration: test_args.a_tr,
    })

    acceleration_value = QuantityCoordinateVector.from_expr(acceleration_subs)

    assert_equal_vectors(acceleration_value, test_args.a_cor)


def test_translation_law(test_args: Args) -> None:
    acceleration_expr = solve_for_vector(law.law, law.translation_acceleration)

    acceleration_subs = acceleration_expr.subs({
        law.relative_acceleration: test_args.a_rel,
        law.force: test_args.f,
        law.mass: test_args.m,
        law.coriolis_acceleration: test_args.a_cor,
    })

    acceleration_value = QuantityCoordinateVector.from_expr(acceleration_subs)

    assert_equal_vectors(acceleration_value, test_args.a_tr)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(mb, test_args.f, test_args.a_cor, test_args.a_tr)
    with raises(TypeError):
        law.calculate_relative_acceleration(100, test_args.f, test_args.a_cor, test_args.a_tr)


def test_bad_force(test_args: Args) -> None:
    fb_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(test_args.m, fb_vector, test_args.a_cor, test_args.a_tr)

    fb_scalar = Quantity(units.newton)
    with raises(ValueError):
        law.calculate_relative_acceleration(test_args.m, fb_scalar, test_args.a_cor, test_args.a_tr)

    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, 100, test_args.a_cor, test_args.a_tr)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, [100], test_args.a_cor, test_args.a_tr)


def test_bad_acceleration(test_args: Args) -> None:
    ab_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, ab_vector, test_args.a_tr)
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, test_args.a_cor, ab_vector)

    ab_scalar = Quantity(units.meter / units.second**2)
    with raises(ValueError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, ab_scalar, test_args.a_tr)
    with raises(ValueError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, test_args.a_cor, ab_scalar)

    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, 100, test_args.a_tr)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, [100], test_args.a_tr)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, test_args.a_cor, 100)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.m, test_args.f, test_args.a_cor, [100])
