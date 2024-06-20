from collections import namedtuple
from pytest import fixture, raises
from sympy import symbols
from symplyphysics import (
    errors,
    units,
    Quantity,
    QuantityVector,
)
from symplyphysics.core.approx import assert_equal_vectors
from symplyphysics.laws.kinematic.vector import acceleration_relative_to_reference_frame as law

Args = namedtuple("Args", "r0 r1 r2 v0 v a dt t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    r0 = QuantityVector([1 * units.meter, 1 * units.meter, -1 * units.meter])
    r1 = QuantityVector([4 * units.meter, 4 * units.meter, 0 * units.meter])
    r2 = QuantityVector([13 * units.meter, 7 * units.meter, 1 * units.meter])
    v0 = QuantityVector([0, 3 * units.meter / units.second, 1 * units.meter / units.second])
    v = QuantityVector([6 * units.meter / units.second, 3 * units.meter / units.second, 1 * units.meter / units.second])
    a = QuantityVector([6 * units.meter / units.second**2, 0, 0])
    dt = Quantity(1 * units.second)
    t = symbols("t")
    return Args(r0=r0, r1=r1, r2=r2, v0=v0, v=v, a=a, dt=dt, t=t)


def test_acceleration_law(test_args: Args) -> None:
    acceleration_ = law.calculate_relative_acceleration(test_args.r0, test_args.r1, test_args.r2, test_args.dt)
    assert_equal_vectors(acceleration_, test_args.a)


def test_velocity_law(test_args: Args) -> None:
    velocity_vector = law.relative_velocity_law(
        test_args.v0.to_base_vector(),
        test_args.a.to_base_vector(),
        test_args.t,
    )
    velocity_ = QuantityVector.from_base_vector(
        velocity_vector,
        subs={test_args.t: test_args.dt},
    )
    assert_equal_vectors(velocity_, test_args.v)


def test_position_law(test_args: Args) -> None:
    position_vector = law.relative_position_law(
        test_args.r0.to_base_vector(),
        test_args.v0.to_base_vector(),
        test_args.a.to_base_vector(),
        test_args.t,
    )
    position_current_ = QuantityVector.from_base_vector(
        position_vector,
        subs={test_args.t: test_args.dt},
    )
    assert_equal_vectors(position_current_, test_args.r1)
    position_after_ = QuantityVector.from_base_vector(
        position_vector,
        subs={test_args.t: 2 * test_args.dt},
    )
    assert_equal_vectors(position_after_, test_args.r2)


def test_bad_positions(test_args: Args) -> None:
    rb_vector = QuantityVector([1 * units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(rb_vector, test_args.r1, test_args.r2, test_args.dt)
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(test_args.r0, rb_vector, test_args.r2, test_args.dt)
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(test_args.r0, test_args.r1, rb_vector, test_args.dt)

    rb_scalar = Quantity(1 * units.meter)
    with raises(AttributeError):
        law.calculate_relative_acceleration(rb_scalar, test_args.r1, test_args.r2, test_args.dt)
    with raises(AttributeError):
        law.calculate_relative_acceleration(test_args.r0, rb_scalar, test_args.r2, test_args.dt)
    with raises(AttributeError):
        law.calculate_relative_acceleration(test_args.r0, test_args.r1, rb_scalar, test_args.dt)

    with raises(TypeError):
        law.calculate_relative_acceleration(100, test_args.r1, test_args.r2, test_args.dt)
    with raises(TypeError):
        law.calculate_relative_acceleration([100], test_args.r1, test_args.r2, test_args.dt)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.r0, 100, test_args.r2, test_args.dt)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.r0, [100], test_args.r2, test_args.dt)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.r0, test_args.r1, 100, test_args.dt)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.r0, test_args.r1, [100], test_args.dt)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_relative_acceleration(test_args.r0, test_args.r1, test_args.r2, tb)
    with raises(TypeError):
        law.calculate_relative_acceleration(test_args.r0, test_args.r1, test_args.r2, 100)
