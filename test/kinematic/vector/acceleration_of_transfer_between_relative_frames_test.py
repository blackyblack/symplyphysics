from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    QuantityVector,
)
from symplyphysics.core.approx import assert_equal_vectors
from symplyphysics.laws.kinematic.vector import acceleration_of_transfer_between_relative_frames as law

Args = namedtuple("Args", "a w0 w1 r t0 t1 t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = QuantityVector([1 * units.meter / units.second**2, 0, 0])
    w0 = QuantityVector([0, 3 * units.radian / units.second, 0])
    w1 = QuantityVector([0, 3.1 * units.radian / units.second, -0.1 * units.radian / units.second])
    r = QuantityVector([3 * units.meter, -1 * units.meter, 2 * units.meter])
    t0 = Quantity(0)
    t1 = Quantity(1 * units.second)
    t = Quantity(0.5 * units.second)
    return Args(a=a, w0=w0, w1=w1, r=r, t0=t0, t1=t1, t=t)


def test_law(test_args: Args) -> None:
    result = law.calculate_transfer_acceleration(test_args.a, test_args.w0, test_args.w1, test_args.r, test_args.t0, test_args.t1, test_args.t)
    a_unit = units.meter / units.second**2
    assert_equal_vectors(result, QuantityVector([-27.0 * a_unit, -0.6 * a_unit, -19.0 * a_unit]), tolerance=2e-2)


def test_bad_acceleration(test_args: Args) -> None:
    ab_vector = QuantityVector([1 * units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_transfer_acceleration(ab_vector, test_args.w0, test_args.w1, test_args.r, test_args.t0, test_args.t1, test_args.t)
    
    ab_scalar = Quantity(1 * units.meter / units.second**2)
    with raises(AttributeError):
        law.calculate_transfer_acceleration(ab_scalar, test_args.w0, test_args.w1, test_args.r, test_args.t0, test_args.t1, test_args.t)
    
    with raises(TypeError):
        law.calculate_transfer_acceleration(100, test_args.w0, test_args.w1, test_args.r, test_args.t0, test_args.t1, test_args.t)
    with raises(TypeError):
        law.calculate_transfer_acceleration([100], test_args.w0, test_args.w1, test_args.r, test_args.t0, test_args.t1, test_args.t)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityVector([1 * units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_transfer_acceleration(test_args.a, wb_vector, test_args.w1, test_args.r, test_args.t0, test_args.t1, test_args.t)
    with raises(errors.UnitsError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, wb_vector, test_args.r, test_args.t0, test_args.t1, test_args.t)

    wb_scalar = Quantity(1 * units.radian / units.second)
    with raises(AttributeError):
        law.calculate_transfer_acceleration(test_args.a, wb_scalar, test_args.w1, test_args.r, test_args.t0, test_args.t1, test_args.t)
    with raises(AttributeError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, wb_scalar, test_args.r, test_args.t0, test_args.t1, test_args.t)

    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a, 100, test_args.w1, test_args.r, test_args.t0, test_args.t1, test_args.t)
    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a, [100], test_args.w1, test_args.r, test_args.t0, test_args.t1, test_args.t)
    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, 100, test_args.r, test_args.t0, test_args.t1, test_args.t)
    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, [100], test_args.r, test_args.t0, test_args.t1, test_args.t)


def test_bad_position_vector(test_args: Args) -> None:
    rb_vector = QuantityVector([1 * units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, test_args.w1, rb_vector, test_args.t0, test_args.t1, test_args.t)

    rb_scalar = Quantity(1 * units.meter)
    with raises(AttributeError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, test_args.w1, rb_scalar, test_args.t0, test_args.t1, test_args.t)

    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, test_args.w1, 100, test_args.t0, test_args.t1, test_args.t)
    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, test_args.w1, [100], test_args.t0, test_args.t1, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, test_args.w1, test_args.r, tb, test_args.t1, test_args.t)
    with raises(errors.UnitsError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, test_args.w1, test_args.r, test_args.t0, tb, test_args.t)
    with raises(errors.UnitsError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, test_args.w1, test_args.r, test_args.t0, test_args.t1, tb)

    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, test_args.w1, test_args.r, 100, test_args.t1, test_args.t)
    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, test_args.w1, test_args.r, test_args.t0, 100, test_args.t)
    with raises(TypeError):
        law.calculate_transfer_acceleration(test_args.a, test_args.w0, test_args.w1, test_args.r, test_args.t0, test_args.t1, 100)
