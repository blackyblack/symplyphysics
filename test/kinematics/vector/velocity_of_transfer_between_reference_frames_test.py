from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.kinematics.vector import velocity_of_transfer_between_reference_frames as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors
from symplyphysics.core.experimental.solvers import solve_for_vector

Args = namedtuple("Args", "vtr v0 w r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    velocity_unit = units.meter / units.second
    vtr = QuantityCoordinateVector([
        3 * velocity_unit,
        -3 * velocity_unit,
        6 * velocity_unit,
    ], CARTESIAN)
    v0 = QuantityCoordinateVector([
        1 * velocity_unit,
        -1 * velocity_unit,
        4 * velocity_unit,
    ], CARTESIAN)
    w = QuantityCoordinateVector([
        Quantity(1 * units.radian / units.second),
        Quantity(0),
        Quantity(-1 * units.radian / units.second),
    ], CARTESIAN)
    r = QuantityCoordinateVector([
        3 * units.meter,
        2 * units.meter,
        -1 * units.meter,
    ], CARTESIAN)
    return Args(vtr=vtr, v0=v0, w=w, r=r)


def test_transfer_law(test_args: Args) -> None:
    result = law.calculate_transfer_velocity(test_args.v0, test_args.w, test_args.r)
    assert_equal_vectors(result, test_args.vtr)


def test_origin_law(test_args: Args) -> None:
    result = solve_for_vector(law.law, law.moving_frame_velocity).subs({
        law.transfer_velocity: test_args.vtr,
        law.angular_velocity: test_args.w,
        law.position_vector: test_args.r,
    })

    result = QuantityCoordinateVector.from_expr(result)

    assert_equal_vectors(result, test_args.v0)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityCoordinateVector([Quantity(1 * units.coulomb), 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_transfer_velocity(vb_vector, test_args.w, test_args.r)

    vb_scalar = Quantity(1 * units.meter / units.second)
    with raises(ValueError):
        law.calculate_transfer_velocity(vb_scalar, test_args.w, test_args.r)

    with raises(TypeError):
        law.calculate_transfer_velocity(100, test_args.w, test_args.r)
    with raises(TypeError):
        law.calculate_transfer_velocity([100], test_args.w, test_args.r)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityCoordinateVector([Quantity(1 * units.coulomb), 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_transfer_velocity(test_args.v0, wb_vector, test_args.r)

    wb_scalar = Quantity(1 * units.radian / units.second)
    with raises(ValueError):
        law.calculate_transfer_velocity(test_args.v0, wb_scalar, test_args.r)

    with raises(TypeError):
        law.calculate_transfer_velocity(test_args.v0, 100, test_args.r)
    with raises(TypeError):
        law.calculate_transfer_velocity(test_args.v0, [100], test_args.r)


def test_bad_position(test_args: Args) -> None:
    rb_vector = QuantityCoordinateVector([Quantity(1 * units.coulomb), 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_transfer_velocity(test_args.v0, test_args.w, rb_vector)

    rb_scalar = Quantity(1 * units.meter)
    with raises(ValueError):
        law.calculate_transfer_velocity(test_args.v0, test_args.w, rb_scalar)

    with raises(TypeError):
        law.calculate_transfer_velocity(test_args.v0, test_args.w, 100)
    with raises(TypeError):
        law.calculate_transfer_velocity(test_args.v0, test_args.w, [100])
