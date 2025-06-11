from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.kinematics.vector import centripetal_acceleration_via_cross_product as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "w r a")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w_unit = units.radian / units.second
    w = QuantityCoordinateVector([0, 1 * w_unit, -1 * w_unit], CARTESIAN)

    r = QuantityCoordinateVector([
        1 * units.meter,
        2 * units.meter,
        -1 * units.meter,
    ], CARTESIAN)

    a_unit = units.meter / units.second**2
    a = QuantityCoordinateVector([-2 * a_unit, -1 * a_unit, -1 * a_unit], CARTESIAN)

    return Args(w=w, r=r, a=a)


def test_cross_law(test_args: Args) -> None:
    result = law.calculate_centripetal_acceleration(test_args.w, test_args.r)
    assert_equal_vectors(result, test_args.a)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityCoordinateVector([1 * units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_centripetal_acceleration(wb_vector, test_args.r)

    wb_scalar = Quantity(1 * units.radian / units.second)
    with raises(ValueError):
        law.calculate_centripetal_acceleration(wb_scalar, test_args.r)

    with raises(TypeError):
        law.calculate_centripetal_acceleration(100, test_args.r)
    with raises(TypeError):
        law.calculate_centripetal_acceleration([100], test_args.r)


def test_bad_position(test_args: Args) -> None:
    rb_vector = QuantityCoordinateVector([1 * units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_centripetal_acceleration(test_args.w, rb_vector)

    rb_scalar = Quantity(1 * units.meter)
    with raises(ValueError):
        law.calculate_centripetal_acceleration(test_args.w, rb_scalar)

    with raises(TypeError):
        law.calculate_centripetal_acceleration(test_args.w, 100)
    with raises(TypeError):
        law.calculate_centripetal_acceleration(test_args.w, [100])
