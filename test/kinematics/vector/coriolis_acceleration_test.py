from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units
from symplyphysics.laws.kinematics.vector import coriolis_acceleration as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "w v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w_unit = units.radian / units.second
    w = QuantityCoordinateVector([1 * w_unit, -1 * w_unit, 0], CARTESIAN)

    v_unit = units.meter / units.second
    v = QuantityCoordinateVector([1 * v_unit, 1 * v_unit, 1 * v_unit], CARTESIAN)
    return Args(w=w, v=v)


def test_law(test_args: Args) -> None:
    result = law.calculate_coriolis_acceleration(test_args.w, test_args.v)

    a_unit = units.meter / units.second**2
    expected = QuantityCoordinateVector([2 * a_unit, 2 * a_unit, -4 * a_unit], CARTESIAN)

    assert_equal_vectors(result, expected)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityCoordinateVector([1 * units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_coriolis_acceleration(wb_vector, test_args.v)

    wb_scalar = 3 * units.radian / units.second
    with raises(ValueError):
        law.calculate_coriolis_acceleration(wb_scalar, test_args.v)

    with raises(TypeError):
        law.calculate_coriolis_acceleration(100, test_args.v)
    with raises(TypeError):
        law.calculate_coriolis_acceleration([100], test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityCoordinateVector([1 * units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_coriolis_acceleration(test_args.w, vb_vector)

    vb_scalar = 3 * units.meter / units.second
    with raises(ValueError):
        law.calculate_coriolis_acceleration(test_args.w, vb_scalar)

    with raises(TypeError):
        law.calculate_coriolis_acceleration(test_args.w, 100)
    with raises(TypeError):
        law.calculate_coriolis_acceleration(test_args.w, [100])
