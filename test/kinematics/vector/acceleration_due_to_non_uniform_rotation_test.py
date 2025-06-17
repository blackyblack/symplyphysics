from collections import namedtuple
from pytest import fixture, raises
from sympy import cos, sin
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.kinematics.vector import acceleration_due_to_non_uniform_rotation as law

from symplyphysics.core.experimental.coordinate_systems import (QuantityCoordinateVector, CARTESIAN,
    CoordinateVector)
from symplyphysics.core.experimental.approx import assert_equal_vectors
from symplyphysics.core.experimental.solvers import vector_equals

Args = namedtuple("Args", "dw dt r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dw = QuantityCoordinateVector([
        0,
        0.001 * units.radian / units.second,
        -0.002 * units.radian / units.second,
    ], CARTESIAN)

    dt = Quantity(0.001 * units.second)

    r = QuantityCoordinateVector([
        1 * units.meter,
        -3 * units.meter,
        0.5 * units.meter,
    ], CARTESIAN)

    return Args(dw=dw, dt=dt, r=r)


def test_law(test_args: Args) -> None:
    result = law.calculate_non_uniform_rotation_acceleration(test_args.dw, test_args.dt,
        test_args.r)

    a_unit = units.meter / units.second**2
    expected = QuantityCoordinateVector([
        -5.5 * a_unit,
        -2.0 * a_unit,
        -1.0 * a_unit,
    ], CARTESIAN)

    assert_equal_vectors(result, expected)


def test_function_law() -> None:
    t = law.time
    w = CoordinateVector([cos(t), 3 * sin(t), 1], CARTESIAN)
    r = CoordinateVector([1, -3, 2], CARTESIAN)

    a_result = law.law.rhs.subs({
        law.angular_velocity(law.time): w,
        law.position_vector: r,
    }).doit()

    a_correct = CoordinateVector([6 * cos(t), 2 * sin(t), 3 * (sin(t) - cos(t))], CARTESIAN)

    assert vector_equals(a_result, a_correct)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_non_uniform_rotation_acceleration(wb_vector, test_args.dt, test_args.r)

    wb_scalar = units.radian / units.second
    with raises(ValueError):
        law.calculate_non_uniform_rotation_acceleration(wb_scalar, test_args.dt, test_args.r)

    with raises(TypeError):
        law.calculate_non_uniform_rotation_acceleration(100, test_args.dt, test_args.r)
    with raises(TypeError):
        law.calculate_non_uniform_rotation_acceleration([100], test_args.dt, test_args.r)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, tb, test_args.r)
    with raises(TypeError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, 100, test_args.r)


def test_bad_position_vector(test_args: Args) -> None:
    rb_vector = QuantityCoordinateVector([units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, test_args.dt, rb_vector)

    rb_scalar = units.meter
    with raises(ValueError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, test_args.dt, rb_scalar)

    with raises(ValueError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, test_args.dt, 100)
    with raises(ValueError):
        law.calculate_non_uniform_rotation_acceleration(test_args.dw, test_args.dt, [100])
