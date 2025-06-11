from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.kinematics.vector import velocity_relative_to_reference_frame as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "r1 r2 dt v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    r1 = QuantityCoordinateVector([
        0,
        units.meter,
        units.meter,
    ], CARTESIAN)

    r2 = QuantityCoordinateVector([
        0.1 * units.meter,
        0.9 * units.meter,
        1.2 * units.meter,
    ], CARTESIAN)

    dt = Quantity(0.1 * units.second)

    v = QuantityCoordinateVector([
        1 * units.meter / units.second,
        -1 * units.meter / units.second,
        2 * units.meter / units.second,
    ], CARTESIAN)

    return Args(r1=r1, r2=r2, dt=dt, v=v)


def test_velocity_law(test_args: Args) -> None:
    result = law.calculate_relative_velocity(test_args.r1, test_args.r2, test_args.dt)
    assert_equal_vectors(result, test_args.v)


# TODO: implement `VectorIntegrate` and add a test case for integrated law


def test_bad_position(test_args: Args) -> None:
    rb_vector = QuantityCoordinateVector([Quantity(1 * units.coulomb), 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_relative_velocity(rb_vector, test_args.r2, test_args.dt)
    with raises(errors.UnitsError):
        law.calculate_relative_velocity(test_args.r1, rb_vector, test_args.dt)

    rb_scalar = Quantity(1 * units.meter)
    with raises(ValueError):
        law.calculate_relative_velocity(rb_scalar, test_args.r2, test_args.dt)
    with raises(ValueError):
        law.calculate_relative_velocity(test_args.r1, rb_scalar, test_args.dt)

    with raises(TypeError):
        law.calculate_relative_velocity(100, test_args.r2, test_args.dt)
    with raises(TypeError):
        law.calculate_relative_velocity(test_args.r1, 100, test_args.dt)
    with raises(TypeError):
        law.calculate_relative_velocity([100], test_args.r2, test_args.dt)
    with raises(TypeError):
        law.calculate_relative_velocity(test_args.r1, [100], test_args.dt)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_relative_velocity(test_args.r1, test_args.r2, tb)
    with raises(TypeError):
        law.calculate_relative_velocity(test_args.r1, test_args.r2, 100)
