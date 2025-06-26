from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, errors, Quantity
from symplyphysics.laws.electricity.vector import electrostatic_force_between_two_charges as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "q1 q2 r f")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q1 = Quantity(1e-10 * units.coulomb)

    q2 = Quantity(-4e-10 * units.coulomb)

    r = QuantityCoordinateVector([
        0.4 * units.milli * units.meter,
        -0.3 * units.milli * units.meter,
        0,
    ], CARTESIAN)

    f = QuantityCoordinateVector([
        -1.15e-3 * units.newton,
        8.63e-4 * units.newton,
        0,
    ], CARTESIAN)

    return Args(q1=q1, q2=q2, r=r, f=f)


def test_force_law(test_args: Args) -> None:
    result = law.calculate_force(test_args.q1, test_args.q2, test_args.r)

    assert_equal_vectors(result, test_args.f)


def test_position_vector_law(test_args: Args) -> None:
    result = law.calculate_position_vector(test_args.q1, test_args.q2, test_args.f)

    assert_equal_vectors(result, test_args.r)


def test_bad_charge(test_args: Args) -> None:
    qb = Quantity(units.candela)

    with raises(errors.UnitsError):
        law.calculate_force(qb, test_args.q2, test_args.r)
    with raises(TypeError):
        law.calculate_force(100, test_args.q2, test_args.r)
    with raises(errors.UnitsError):
        law.calculate_force(test_args.q1, qb, test_args.r)
    with raises(TypeError):
        law.calculate_force(test_args.q1, 100, test_args.r)

    with raises(errors.UnitsError):
        law.calculate_position_vector(qb, test_args.q2, test_args.f)
    with raises(TypeError):
        law.calculate_position_vector(100, test_args.q2, test_args.f)
    with raises(errors.UnitsError):
        law.calculate_position_vector(test_args.q1, qb, test_args.f)
    with raises(TypeError):
        law.calculate_position_vector(test_args.q1, 100, test_args.f)


def test_bad_position_vector(test_args: Args) -> None:
    rb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_force(test_args.q1, test_args.q2, rb_vector)

    rb_scalar = Quantity(units.meter)
    with raises(ValueError):
        law.calculate_force(test_args.q1, test_args.q2, rb_scalar)

    with raises(TypeError):
        law.calculate_force(test_args.q1, test_args.q2, 100)
    with raises(TypeError):
        law.calculate_force(test_args.q1, test_args.q2, [100])


def test_bad_force(test_args: Args) -> None:
    fb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_position_vector(test_args.q1, test_args.q2, fb_vector)

    fb_scalar = units.newton
    with raises(ValueError):
        law.calculate_position_vector(test_args.q1, test_args.q2, fb_scalar)

    with raises(TypeError):
        law.calculate_position_vector(test_args.q1, test_args.q2, 100)
    with raises(TypeError):
        law.calculate_position_vector(test_args.q1, test_args.q2, [100])
