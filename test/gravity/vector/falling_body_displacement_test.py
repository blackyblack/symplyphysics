from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.gravity.vector import falling_body_displacement as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "t v0 w g")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(5 * units.second)

    v_unit = units.meter / units.second
    v0 = QuantityCoordinateVector([v_unit, v_unit, v_unit], CARTESIAN)

    w = QuantityCoordinateVector([0, -7.27e-5 * units.radian / units.second, 0], CARTESIAN)

    g = QuantityCoordinateVector([0, 0, -9.81 * units.meter / units.second**2], CARTESIAN)

    return Args(t=t, v0=v0, w=w, g=g)


def test_law(test_args: Args) -> None:
    result = law.calculate_displacement(test_args.t, test_args.v0, test_args.w, test_args.g)

    expected = QuantityCoordinateVector([
        4.97 * units.meter,
        5.00 * units.meter,
        -118.0 * units.meter,
    ], CARTESIAN)

    assert_equal_vectors(result, expected, relative_tolerance=4e-3)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_displacement(tb, test_args.v0, test_args.w, test_args.g)
    with raises(TypeError):
        law.calculate_displacement(100, test_args.v0, test_args.w, test_args.g)


def test_bad_linear_velocity(test_args: Args) -> None:
    vb_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_displacement(test_args.t, vb_vector, test_args.w, test_args.g)

    vb_scalar = Quantity(units.meter / units.second)
    with raises(ValueError):
        law.calculate_displacement(test_args.t, vb_scalar, test_args.w, test_args.g)

    with raises(TypeError):
        law.calculate_displacement(test_args.t, 100, test_args.w, test_args.g)
    with raises(TypeError):
        law.calculate_displacement(test_args.t, [100], test_args.w, test_args.g)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_displacement(test_args.t, test_args.v0, wb_vector, test_args.g)

    wb_scalar = Quantity(units.radian / units.second)
    with raises(ValueError):
        law.calculate_displacement(test_args.t, test_args.v0, wb_scalar, test_args.g)

    with raises(TypeError):
        law.calculate_displacement(test_args.t, test_args.v0, 100, test_args.g)
    with raises(TypeError):
        law.calculate_displacement(test_args.t, test_args.v0, [100], test_args.g)


def test_bad_acceleration(test_args: Args) -> None:
    gb_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_displacement(test_args.t, test_args.v0, test_args.w, gb_vector)

    gb_scalar = Quantity(units.meter / units.second**2)
    with raises(ValueError):
        law.calculate_displacement(test_args.t, test_args.v0, test_args.w, gb_scalar)

    with raises(TypeError):
        law.calculate_displacement(test_args.t, test_args.v0, test_args.w, 100)
    with raises(TypeError):
        law.calculate_displacement(test_args.t, test_args.v0, test_args.w, [100])
